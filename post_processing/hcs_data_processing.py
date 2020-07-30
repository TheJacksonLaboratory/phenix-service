import os
import sys
import typing
import logging
import argparse
import numpy as np
import xlrd
import pandas as pd
from pathlib import Path

import utils
from plate_plotting import plate_qc, plot_randomization, plot_measurements, PdfPages


class HighContentScreen:
    def __init__(self, input_form_path, overwrite=False):
        self._allowed_measurements = {
            "spectramax": self._load_spectramax,
            "phenix": self._load_phenix,
        }
        self.has_randomization = False
        self.input_form_path = input_form_path
        self.overwrite_allowed = overwrite
        self.parse_input_form()
        self.data = []
        logger.debug("HCS object initialized successfully.")

    @staticmethod
    def validate_hcs_excel_file(excel_file):
        log_assert(
            len(excel_file.sheet_names) == 2, "Excel file malformed.  Expected 2 sheets"
        )

    def _parse_screen_metadata(self, excel_file):
        raw = excel_file.parse("Description", index_col=0, header=None).fillna("")
        description = raw.iloc[:18, :1].squeeze()
        new_index = description.index.dropna()
        description = description.loc[new_index].to_dict()

        dispensing = raw.iloc[19:22]
        dispensing = dispensing.rename(columns=dispensing.iloc[0]).drop(
            dispensing.index[0]
        )
        dilution_points = raw.iloc[22, 0]
        dilution_constant = raw.iloc[23, 0]
        notes = raw.iloc[24, 0]

        description["Serial dilution constant"] = dilution_constant
        description["Serial dilutions"] = dilution_points
        description["Notes"] = notes

        self.metadata = description
        self.serial_dilution_constant = dilution_constant
        self.serial_dilution_series_length = dilution_points
        self.dispensing = dispensing

    def _parse_screen_plate_variables(self, excel_file):
        variables = excel_file.parse("Imaged Plates", index_col=0, header=0)
        variable_names = variables["Variable Name"]
        valid_names = ~variable_names.isnull()
        n_variables = np.sum(valid_names)
        if n_variables == 0:
            logger.info("No daughter plate info defined.")
            variables = None
        else:
            variables = variables.loc[valid_names]
            variables.set_index(variables.columns[0], drop=True, inplace=True)
            logger.debug(f"Found the following plate variables defined:\n{variables}")

        self.variables = variables

    def parse_input_form(self):
        logger.debug("Parsing HCS excel file...")
        input_excel = pd.ExcelFile(self.input_form_path)
        self.validate_hcs_excel_file(input_excel)
        self._parse_screen_metadata(input_excel)
        self._parse_screen_plate_variables(input_excel)
        logger.debug("HCS excel file read successfully.")

    def construct_dilution_series(self):
        dilution_series = (1 / self.serial_dilution_constant) ** np.arange(
            self.serial_dilution_series_length
        )

        dilutions = self.dispensing.iloc[1, :].values[:, None] * (dilution_series)
        dilution_plate = np.zeros((8, 12))
        dilution_plate[:, : self.serial_dilution_series_length] = dilutions
        dilution_data = utils.construct_96(
            utils.flatten_plate_map(dilution_plate, colwise=False),
            "Concentration",
            colwise=False,
        )
        self.dilution_series = dilution_data
        logger.debug("Created dilution series data sucessfully.")

    def construct_drug_series(self):
        self.drug_series = utils.construct_96(
            utils.flatten_plate_map(
                np.repeat(self.dispensing.iloc[0, :].values, 12), colwise=False
            ),
            "Drug",
            colwise=False,
        )
        logger.debug("Created drug name data sucessfully.")

    def construct_randomization(self):
        sheet_names = pd.ExcelFile(self.randomization_file).sheet_names
        quads = []
        for k, name in enumerate(sheet_names):
            with open(os.devnull, "w") as devnull:
                wb = xlrd.open_workbook(self.randomization_file, logfile=devnull)
                quad_data = (
                    pd.read_excel(wb, sheet_name=name, index_col=0, engine="xlrd")
                    .dropna()
                    .iloc[:, 0]
                    .values
                )
            quad = utils.convert_96_to_384(
                quad_data, name="Source well 96", quad=k, colwise=True
            )
            quads.append(quad)

        randomization_data = pd.concat(quads, axis=0)
        (
            randomization_data["Source row 96"],
            randomization_data["Source col 96"],
        ) = utils.split_row_col(randomization_data["Source well 96"]).values.T
        (
            randomization_data["Row 384"],
            randomization_data["Col 384"],
        ) = utils.split_row_col(randomization_data.index).values.T
        self.data.append(randomization_data)
        self.has_randomization = True
        self.randomization_mapping = randomization_data.iloc[:, 0].to_dict()
        logger.debug("Randomization mapping created successfully.")

    def map_randomization(self):
        if not self.has_randomization:
            logger_exception(
                "Must load randomization before trying to map randomization."
            )
            exit(2)

        dilution_series_384 = self.dilution_series.reindex(
            self.randomization_mapping.values()
        )
        dilution_series_384.index = self.randomization_mapping.keys()
        drug_series_384 = self.drug_series.reindex(self.randomization_mapping.values())
        drug_series_384.index = self.randomization_mapping.keys()
        self.data.append(dilution_series_384)
        self.data.append(drug_series_384)
        logger.debug("Dilution and drug data mapped using randomization mapping.")

    @staticmethod
    def _load_spectramax(spectramax_file, custom_name, ignored):
        data = (
            pd.read_table(
                spectramax_file,
                encoding="utf-16",
                engine="python",
                sep="\t",
                skiprows=2,
                skipfooter=2,
            )
            .iloc[:16, 2:26]
            .values
        )
        flat_data = utils.flatten_plate_map(data, colwise=False)
        logger.debug(f"[{custom_name}] data constructed successfully.")
        data = utils.construct_384(flat_data, custom_name, colwise=False)
        return data, None

    @staticmethod
    def _load_phenix(phenix_file, custom_name, phenix_columns):
        data = pd.read_table(phenix_file, index_col=2)
        cols = utils.parse_column_spec(phenix_columns)
        logger.debug(f"Will use [{custom_name}] data columns: [{cols}]")
        try:
            data = data.loc[:, cols]
        except:
            data = data.iloc[:, cols]
        logger.debug("Phenix data constructed successfully.")
        if len(data.shape) == 1:
            data = data.to_frame()
        data.columns = f"{custom_name} - " + data.columns
        return data, data.columns.tolist()

    def load_measurements(self):
        measurements = []
        self.measurements = []
        for m_name, (m_file, m_cols) in self.measurement_files.items():
            assert_path_exists(m_name, m_file)
            is_valid = False
            m_key = None
            for key_name in self._allowed_measurements.keys():
                if key_name in m_name.lower():
                    is_valid = True
                    m_key = key_name
            if not is_valid:
                logger.warn(f"Measurement {m_name} not configured. Skipping.")
            logger.debug(f"Trying to load {m_name} from file [{m_file}]")
            m_data, new_names = self._allowed_measurements[m_key](
                m_file, m_name, m_cols
            )
            measurements.append(m_data)
            if new_names:
                self.measurements += new_names
            else:
                self.measurements.append(m_name)

        logger.debug("All measurements loaded successfully.")
        self.data.append(pd.concat(measurements, axis=1))

    def register_data(self, measurements={}, randomization=None):
        if measurements:
            self.measurement_files = measurements

        if randomization:
            self.randomization_file = randomization

    def aggregate_data(self):
        merged = pd.concat(self.data, axis=1)
        merged.index.name = "Well 384"
        merged = utils.sort_index(merged)

        for col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="ignore")
        for col in merged.columns[merged.dtypes == "object"]:
            merged[col] = merged[col].astype("category")

        logger.debug("Master data series constructed successfully.")
        return merged

    def incorporate_plate_variables(self, merged):
        if self.variables is not None:
            n_plates = self.variables.shape[1]
            expanded = pd.concat([merged] * n_plates, axis=0)
            repeat = lambda x: np.repeat(x, merged.shape[0])
            expanded["plate"] = repeat(np.arange(n_plates, dtype=int) + 1)
            for key, row in self.variables.iterrows():
                expanded[key] = repeat(row.values)
            logger.debug("Expanded data series constructed successfully.")
        else:
            expanded = merged.copy()
            logger.debug("No plate variables, so not expanding data series.")
        return expanded

    def pipeline(self):
        self.construct_dilution_series()
        self.construct_drug_series()
        self.construct_randomization()

        self.map_randomization()
        self.load_measurements()
        aggregated = self.aggregate_data()
        finalized = self.incorporate_plate_variables(aggregated)
        self.output = finalized
        return finalized

    def save_output(self, output_file):
        if not self.overwrite_allowed:
            assert_path_does_not_exist("Output file", output_file)
        self.output.to_csv(output_file)
        logger.debug(f"Output saved to [{output_file}].")

    def plot_plate_maps(self, plot_file):
        logger.debug(f"Attempting to plot plate data...")
        if not self.overwrite_allowed:
            assert_path_does_not_exist("Plot file", plot_file)
        figs = [
            plate_qc(self.output),
            plot_randomization(self.output),
            plot_measurements(self.output, measurement_cols=self.measurements),
        ]
        logger.debug(f"Plots generated")
        pdf = PdfPages(plot_file)
        for fig in figs:
            pdf.savefig(fig)
        pdf.close()
        logger.debug(f"Plots saved to [{plot_file}].")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse related files for use with the Single Cell Biology "
            "Lab Opera Phenix High Content Screening platform"
        ),
        prog=__file__,
    )

    parser.add_argument(
        "-i",
        dest="hcs_input_file",
        type=Path,
        required=True,
        help="Path to HCS input form",
    )

    parser.add_argument(
        "-r",
        dest="randomization_file",
        type=Path,
        required=True,
        help="Path to randomization CSV file",
    )

    parser.add_argument(
        "-m",
        "--measurement",
        dest="measurements",
        nargs="+",
        action="append",
        required=False,
        help=(
            "Measurement specification: measurement_name measurement_file "
            "[columns_to_use]. If [columns_to_use] is not given, will use "
            " all columns presented."
        ),
    )

    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        required=True,
        help="Path to save output dataframe",
    )

    parser.add_argument(
        "-p",
        dest="plot_pdf_file",
        type=Path,
        default=None,
        help="Path to save plots [None = default]",
    )

    parser.add_argument(
        "-f", "--force", action="store_true", help="Overwite output files"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print extra information"
    )

    args = parser.parse_args()
    if args.measurements:
        args.measurements = dict(
            (item[0], [Path(item[1]), None if len(item) == 2 else item[2]])
            for item in args.measurements
        )

    return args


def log_assert(assertion, message):
    try:
        assert assertion
    except AssertionError as err:
        logger.exception(message)
        raise err


def assert_path_exists(name: str, path: Path) -> None:
    if path is not None and not path.exists():
        logger.error(f"{name.replace('_',' ').capitalize()}: [{path}] does not exist!")
        print(path)
        print(path.exists())
        exit(2)


def assert_path_does_not_exist(name: str, path: Path) -> None:
    if path is not None and path.exists():
        logger.error(f"{name.replace('_',' ').capitalize()}: [{path}] already exist!")
        exit(2)


if __name__ == "__main__":
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format="%(asctime)s - %(name)s: %(levelname)s:  %(message)s",
    )
    logger = logging.getLogger(__file__)

    args = parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    logger.debug(f"Parsed arguments: {vars(args)}")

    hcs = HighContentScreen(args.hcs_input_file, overwrite=args.force)
    hcs.register_data(
        measurements=args.measurements, randomization=args.randomization_file,
    )
    aggregated = hcs.pipeline()
    hcs.save_output(args.output_file)

    if args.plot_pdf_file:
        hcs.plot_plate_maps(args.plot_pdf_file)
