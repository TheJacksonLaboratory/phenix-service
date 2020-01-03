import sys
import json
import typing
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


plate_params = {
    96: {"width": 12, "height": 8},  # +4 due to name + space + top/bottom indexes
    384: {"width": 24, "height": 16},
}


def log_assert(assertion, message):
    try:
        assert assertion
    except AssertionError as err:
        logger.exception(message)
        raise err


def parse_description(excel_file: pd.ExcelFile) -> dict:
    raw_description = excel_file.parse("Description", index_col=0)

    description = raw_description.iloc[:, :1].squeeze()
    new_index = description.index.dropna()
    description = description.loc[new_index]
    return description.to_dict()


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


def parse_phenix_metadata(fname: str) -> typing.Tuple[dict, pd.DataFrame]:

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse a 96-well or 384-well plate map for use with the Single Cell Biology "
            "Lab Opera Phenix"
        )
    )

    parser.add_argument("plate_map_path", help="Path to the Excel file")
    parser.add_argument(
        "output_metadata_path", help="Path where to write metadata.json"
    )

    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    json_path = Path(args.output_metadata_path)
    log_assert(
        json_path.parent.exists(), f"Directory containing {json_path} doesn't exist!"
    )

    metadata, plate_dataframe = parse_phenix_metadata(args.plate_map_path)
    metadata.update(plate_dataframe.to_dict(orient="index"))

    with open(str(json_path), "w") as fout:
        json.dump(metadata, fout)

    logger.info(f"Successfully wrote metadat to {json_path}")


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
