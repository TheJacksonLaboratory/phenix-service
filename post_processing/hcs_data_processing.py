import sys
import json
import typing
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse related files for use with the Single Cell Biology "
            "Lab Opera Phenix High Content Screening platform"
        ),
        prog=__file__,
    )
    add_arg = lambda a, **v: parser.add_argument(a, **v)

    add_arg(
        "-i",
        dest="hcs_input_file",
        type=Path,
        required=True,
        help="Path to HCS input form",
    )
    add_arg(
        "-r",
        dest="randomization_file",
        type=Path,
        required=True,
        help="Path to randomization CSV file",
    )
    add_arg(
        "-s",
        dest="spectramax_file",
        type=Path,
        required=True,
        help="Path to Spectra Max export file",
    )
    add_arg(
        "-o",
        dest="output_file",
        type=Path,
        required=True,
        help="Path to Spectra Max export file",
    )
    add_arg(
        "-p",
        dest="plot_pdf_file",
        type=Path,
        default=None,
        help="Path to Spectra Max export file",
    )

    return parser.parse_args()


def log_assert(assertion, message):
    try:
        assert assertion
    except AssertionError as err:
        logger.exception(message)
        raise err


def parse_plates(excel_file: pd.ExcelFile, n_variables: int) -> list:
    plate_size = int(excel_file.sheet_names[-1].split()[0])
    width, height = plate_params[plate_size].values()

    plate_sheet = excel_file.parse(excel_file.sheet_names[-1], header=None)
    plates = [
        plate_sheet.iloc[
            (2 + 4 * k) + k * height : (2 + 4 * k) + (k + 1) * height, 1 : 1 + width + 1
        ]
        for k in range(n_variables)
    ]
    for plate in plates:
        plate.set_index(1, inplace=True, drop=True)
        plate.index.name = None
        plate.columns = range(1, width + 1)
    return plates


def flatten_plate(plate: pd.DataFrame, title: str = None) -> pd.Series:
    flat = plate.stack(dropna=False)
    flat.index = flat.index.map("{0[0]}{0[1]}".format)
    if title:
        flat.name = title
    return flat


def parse_phenix_metadata(path: Path) -> typing.Tuple[dict, pd.DataFrame]:
    valid_suffixes = (".xls", ".xlsx", ".xltx")
    fpath = Path(fname)
    log_assert(fpath.exists(), f"File {fpath} doesn't exist!")
    log_assert(
        fpath.suffix in valid_suffixes,
        f"File {fpath} doesn't have one of the following file extensions "
        f"[{' '.join(valid_suffixes)}]",
    )

    excel_file = pd.ExcelFile(fname)
    log_assert(
        len(excel_file.sheet_names) == 3, "Excel file malformed.  Expected 3 sheets"
    )

    description = parse_description(excel_file)

    _variable_names = excel_file.parse("Variables", index_col=0, header=0)[
        "Variable Name"
    ]
    _valid_names = ~_variable_names.isnull()
    n_variables = np.sum(_valid_names)
    variable_names = _variable_names.loc[_valid_names].values
    if n_variables == 0:
        logger.error("No variables defined in variable sheet!")
        exit()

    plates = parse_plates(excel_file, n_variables)
    flat_plates = [
        flatten_plate(plate, variable_name)
        for plate, variable_name in zip(plates, variable_names)
    ]
    plate_map = pd.concat(flat_plates, axis=1)

    return description, plate_map


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


def parse_description(excel_file: pd.ExcelFile) -> typing.Tuple[dict, pd.DataFrame]:
    raw_description = excel_file.parse("Description", index_col=0, header=None).fillna(
        ""
    )

    description = raw_description.iloc[:18, :1].squeeze()
    new_index = description.index.dropna()
    description = description.loc[new_index].to_dict()

    dispensing = raw_description.iloc[19:22]
    dispensing = dispensing.rename(columns=dispensing.iloc[0]).drop(dispensing.index[0])
    dilution_points = raw_description.iloc[22, 0]
    dilution_constant = raw_description.iloc[23, 0]
    notes = raw_description.iloc[24, 0]

    description["Serial dilution constant"] = dilution_constant
    description["Serial dilutions"] = dilution_points
    description["Notes"] = notes

    return description, dispensing


def parse_daughter_plate_specs(
    excel_file: pd.ExcelFile,
) -> typing.Union[pd.DataFrame, None]:
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

    return variables


def construct_dilution_plate(dispensing, n_dilutions, dilution_constant=10 ** (1 / 2)):
    dilutions = dispensing.iloc[1, :].values[:, None] * (
        (1 / dilution_constant) ** np.arange(n_dilutions)
    )
    dilution_plate = pd.DataFrame(
        0, index=list("ABCDEFGH"), columns=pd.RangeIndex(1, 13)
    )
    dilution_plate.iloc[:, :n_dilutions] = dilutions
    dilution_plate = dilution_plate.stack().to_frame()
    dilution_plate.index = dilution_plate.index.map("{0[0]}{0[1]}".format)
    dilution_plate.columns = ["drug_concentration"]
    dilution_plate["drug"] = np.repeat(dispensing.iloc[0, :].values, 12)
    return dilution_plate


def parse_randomization_file(randomization_file):
    sheet_names = pd.ExcelFile(randomization_file).sheet_names

    sheets = []
    for name in sheet_names:
        sheet = (
            pd.read_excel(randomization_file, sheet_name=name, index_col=0)
            .dropna()
            .iloc[:, 0]
            .values
        )
        sheets.append(sheet.reshape(8, -1))

    data_384 = np.block([sheets[:2], sheets[2:]])
    df_384 = (
        pd.DataFrame(
            data_384, index=list("ABCDEFGHIJKLMNOP"), columns=pd.RangeIndex(1, 25)
        )
        .stack()
        .to_frame()
    )
    df_384.index = df_384.index.map("{0[0]}{0[1]}".format)
    df_384.columns = ["source_well_96"]
    return df_384


def parse_spectramax_file(spectramax_file, encoding="utf-16"):
    """
    This file has a bizarre format, is UTF-16 encoded (why?), and due to being
    generated on Windows has carriage returns as line endings.
    Due to this, the parsing is pretty fragile.
    """
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
        .stack()
        .values
    )
    return data


def construct_dataframe(dilution, randomization_file, spectramax_file):
    df_384 = parse_randomization_file(randomization_file)
    spectra = parse_spectramax_file(spectramax_file)

    dilution.index.name = "source_well_96"
    final = dilution.loc[df_384.source_well_96, :].reset_index()
    final.index = df_384.index
    final.index.name = "well_384"

    final["spectramax"] = parse_spectramax_file(spectramax_file)

    final[["source_row_96", "source_col_96"]] = final.source_well_96.str.extract(
        "([A-Z])(\d+)"
    )
    final[["row_384", "col_384"]] = final.index.to_series().str.extract(
        "([A-Z])(\d+)", expand=True
    )

    for col in final.columns:
        final[col] = pd.to_numeric(final[col], errors="ignore")
    for col in final.columns[final.dtypes == "object"]:
        final[col] = final[col].astype("category")

    final["drug_code"] = final.drug.cat.codes
    final["source_y"] = 9 - final.source_row_96.cat.codes
    final["source_x"] = final.source_col_96
    final["y"] = 16 - final.row_384.cat.codes
    final["x"] = final.col_384

    final["quad"] = \
        1 + 2 * ((final["y"] - 7.5)<0).astype(int) + ((final["x"] - 11.5)>0).astype(int)

    return final


def construct_daughter_dataframe(final, daughter_info):
    n_plates = daughter_info.shape[1]
    expanded = pd.concat([final] * n_plates, axis=0)
    repeat = lambda x: np.repeat(x, final.shape[0])
    expanded["plate"] = repeat(np.arange(n_plates, dtype=int) + 1)
    for key, row in daughter_info.iterrows():
        expanded[key] = repeat(row.values)
    return expanded


def main(args: argparse.Namespace) -> None:
    vargs = vars(args)
    for name in ("output_file", "plot_pdf_file"):
        assert_path_does_not_exist(name, vargs[name])
    for name in ("hcs_input_file", "randomization_file", "spectramax_file"):
        assert_path_exists(name, vargs[name])

    input_excel = pd.ExcelFile(args.hcs_input_file)
    log_assert(
        len(input_excel.sheet_names) == 2, "Excel file malformed.  Expected 2 sheets"
    )
    description, dispensing = parse_description(input_excel)

    daughter_plate_info = parse_daughter_plate_specs(input_excel)

    dilution_plate = construct_dilution_plate(
        dispensing, description["Serial dilutions"], description["Serial dilution constant"]
    )

    dataframe = construct_dataframe(
        dilution_plate, args.randomization_file, args.spectramax_file,
    )

    daughter_dataframe = construct_daughter_dataframe(dataframe, daughter_plate_info)

    if args.plot_pdf_file is not None:
        plot_plates()

    dataframe.to_csv(args.output_file)
    daughter_outfile = (
        args.output_file.parent
        / f"{args.output_file.stem}-expanded{args.output_file.suffix}"
    )
    daughter_dataframe.to_csv(daughter_outfile)


if __name__ == "__main__":
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s: %(funcName)s - %(levelname)s:  %(message)s",
    )
    logger = logging.getLogger(__file__)

    args = parse_args()
    logger.debug(f"Parsed arguments: {vars(args)}")
    main(args)
