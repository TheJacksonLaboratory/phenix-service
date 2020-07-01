import sys
import json
import typing
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import product

from string import ascii_uppercase

COLS_96 = range(1, 13)
ROWS_96 = ascii_uppercase[:8]
COLS_384 = range(1, 25)
ROWS_384 = ascii_uppercase[:16]
INDEX_96_C = pd.Index([f"{r}{c}" for r,c in product(ROWS_96, COLS_96)])
INDEX_96_F = pd.Index(INDEX_96_C.values.reshape(len(ROWS_96), -1).ravel(order="F"))
INDEX_384_C = pd.Index([f"{r}{c}" for r,c in product(ROWS_384, COLS_384)])
INDEX_384_F = pd.Index(INDEX_384_C.values.reshape(len(ROWS_384), -1).ravel(order="F"))


assert len(INDEX_96_C) == len(INDEX_96_F) == 96
assert len(INDEX_384_C) == len(INDEX_384_F) == 384

def well_sort(item):
    return (item[0], int(item[1:]))

def sort_index(series):
    return series.reindex(sorted(series.index, key=well_sort))

class HighContentScreen:

    def __init__(self, input_form_path, overwrite=False):
        self._allowed_measurements = {"spectramax": self._load_spectramax}
        self.has_randomization = False
        self.input_form_path = input_form_path
        self.overwrite_allowed = overwrite
        self.parse_input_form()
        self.data = []
        logger.debug("HCS object initialized successfully.")

    @staticmethod
    def flatten_plate_map(data, colwise=False):
        ravel = "F" if colwise else "C"
        return data.ravel(order=ravel)

    @staticmethod
    def construct_96(data_96, name, colwise=False):
        index = INDEX_96_F if colwise else INDEX_96_C
        return pd.Series(data_96, name=name, index=index)

    @staticmethod
    def construct_384(data_384, name, colwise=False):
        index = INDEX_384_F if colwise else INDEX_384_C
        return pd.Series(data_384, name=name, index=index)

    @staticmethod
    def index_by_quad(quad, colwise=False):
        assert quad < 4
        ravel = "F" if colwise else "C"
        source = INDEX_384_C.values.reshape(len(ROWS_384), -1)
        i, j = quad//2, quad%2
        print(i, j)
        return pd.Index(source[i::2, j::2].ravel(order=ravel))

    @staticmethod
    def convert_96_to_384(data_96, name, quad=None, colwise=False):
        if quad is not None and len(data_96) > 96:
            raise ValueError(f"96 well data for quad={quad} has more than 96 values")
        elif quad is None and len(data_96) == 96:
            raise ValueError(f"96 well data with no quad specified has only 96 values")

        if quad is None:
            index = pd.Index(sum([HighContentScreen.index_by_quad(k, colwise).tolist() for k in range(4)], []))
            quads = np.repeat(range(4), 96)
        else:
            index = HighContentScreen.index_by_quad(quad, colwise)
            quads = np.ones(96, dtype=np.uint8) * (quad + 1)

        s = pd.DataFrame(data_96, index=index, columns=[name])
        s["Quadrant"] = quads

        return s

    @staticmethod
    def split_row_col(s):
        if isinstance(s, list) or isinstance(s, pd.Index):
            s = pd.Series(s)
        return s.str.extract("([A-Z])(\d+)", expand=True)

    @staticmethod
    def validate_hcs_excel_file(excel_file):
        log_assert(
            len(excel_file.sheet_names) == 2, "Excel file malformed.  Expected 2 sheets"
        )

    def _parse_screen_metadata(self, excel_file):
        raw = excel_file.parse("Description", index_col=0, header=None).fillna(
            ""
        )
        description = raw.iloc[:18, :1].squeeze()
        new_index = description.index.dropna()
        description = description.loc[new_index].to_dict()

        dispensing = raw.iloc[19:22]
        dispensing = dispensing.rename(columns=dispensing.iloc[0]).drop(dispensing.index[0])
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
            logger.debug(f"Found the following plate variables defined: {variables}")

        self.variables = variables

    def parse_input_form(self):
        logger.debug("Parsing HCS excel file...")
        input_excel = pd.ExcelFile(self.input_form_path)
        self.validate_hcs_excel_file(input_excel)
        self._parse_screen_metadata(input_excel)
        self._parse_screen_plate_variables(input_excel)
        logger.debug("HCS excel file read successfully.")

    def construct_dilution_series(self):
        dilution_series = (1 / self.serial_dilution_constant) ** np.arange(self.serial_dilution_series_length)

        dilutions = self.dispensing.iloc[1, :].values[:, None] * (
            dilution_series
        )
        dilution_plate = np.zeros((8, 12))
        dilution_plate[:, :self.serial_dilution_series_length] = dilutions
        dilution_data = self.construct_96(
            self.flatten_plate_map(dilution_plate, colwise=False),
            "Concentration", colwise=False
        )
        self.dilution_series = dilution_data
        logger.debug("Created dilution series data sucessfully.")

    def construct_drug_series(self):
        self.drug_series = self.construct_96(
            self.flatten_plate_map(
                np.repeat(self.dispensing.iloc[0,:].values, 12), colwise=False
            ), "Drug", colwise=False
        )
        logger.debug("Created drug name data sucessfully.")

    def construct_randomization(self):
        sheet_names = pd.ExcelFile(self.randomization_file).sheet_names
        quads = []
        for k, name in enumerate(sheet_names):
            quad_data = pd.read_excel(self.randomization_file, sheet_name=name,
                    index_col=0).dropna().iloc[:,0].values
            quad = self.convert_96_to_384(
                quad_data, name="Source well 96", quad=k, colwise=True
            )
            quads.append(quad)

        randomization_data = pd.concat(quads, axis=0)
        randomization_data["Source row 96"], randomization_data["Source col 96"] = self.split_row_col(randomization_data["Source well 96"]).values.T
        randomization_data["Row 384"], randomization_data["Col 384"] = self.split_row_col(randomization_data.index).values.T
        self.data.append(randomization_data)
        self.has_randomization = True
        self.randomization_mapping = randomization_data.iloc[:,0].to_dict()
        logger.debug("Randomization mapping created successfully.")

    def map_randomization(self):
        if not self.has_randomization:
            logger_exception("Must load randomization before trying to map randomization.")
            exit(2)

        dilution_series_384 = self.dilution_series.reindex(self.randomization_mapping.values())
        dilution_series_384.index = self.randomization_mapping.keys()
        drug_series_384 = self.drug_series.reindex(self.randomization_mapping.values())
        drug_series_384.index = self.randomization_mapping.keys()
        self.data.append(dilution_series_384)
        self.data.append(drug_series_384)
        logger.debug("Dilution and drug data mapped using randomization mapping.")

    @staticmethod
    def _load_spectramax(spectramax_file):
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
        flat_data = HighContentScreen.flatten_plate_map(data, colwise=False)
        logger.debug("Spectramax data constructed successfully.")
        return HighContentScreen.construct_384(flat_data, "spectramax", colwise=False)

    def load_measurements(self):
        measurements = []
        for m_name, m_file in self.measurement_files.items():
            assert_path_exists(m_name, m_file)
            if m_name not in self._allowed_measurements.keys():
                logger.warn(f"Measurement {m_name} not configured. Skipping.")
            logger.debug(f"Trying to load {m_name} from file [{m_file}]")
            measurements.append(self._allowed_measurements[m_name](m_file))
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
        merged = sort_index(merged)

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

    def plot_plate_maps(self, output_file):
        raise NotImplementedError("Don't have plotting yet")


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
        "-m", "--measurement",
        dest="measurements",
        nargs=2, action="append",
        required=False,
        help="Path to Spectra Max export file",
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
        "-f", "--force", action="store_true",
        help="Overwite output files"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print extra information"
    )

    args = parser.parse_args()
    if args.measurements:
        args.measurements = dict((item[0], Path(item[1])) for item in args.measurements)

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
        measurements=args.measurements,
        randomization=args.randomization_file
    )
    aggregated = hcs.pipeline()
    hcs.save_output(args.output_file)
