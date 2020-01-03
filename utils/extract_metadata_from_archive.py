#!/usr/bin/env python
"""

"""
import os
import sys
import json
import pathlib
import logging
import argparse

import xml.etree.ElementTree as ET


VALID_DATATYPES = {
    "measurement": {"fields": ["PlateName", "UserName", "MeasurementID",
    "TargetTemperature", "TargetCO2"]},
    "analysissequence": {"fields": []},
    "experiment": {"fields": []},
    "assaylayout": {"fields": []},
}
TAG_PREFIX = "{http://www.perkinelmer.com/PEHH/HarmonyV5}"

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger(__name__)


def search_xml_for_key(xml_path, key):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    found_tag = root.find(f".//{TAG_PREFIX}{key}")
    value = None
    if found_tag is not None:
        value = found_tag.text
    return value

parser = argparse.ArgumentParser()
parser.add_argument("archive_path", type=pathlib.Path)
args = parser.parse_args()

if not args.archive_path.exists():
    logger.error(f"{args.archive_path} doesn't exist!")
    sys.exit(1)

xml_dir = args.archive_path / "XML"
if not xml_dir.exists():
    logger.critical(f"{args.archive_path} is empty! Skipping.")
    sys.exit(1)

archive_metadata = {}
for datatype_dir in xml_dir.iterdir():
    datatype = datatype_dir.name.lower()

    datatype_obj = VALID_DATATYPES.get(datatype, None)
    if not datatype_obj:
        continue

    keys_to_search = datatype_obj["fields"]
    for search_key in keys_to_search:
        for xml_file in datatype_dir.glob("*.xml"):
            if xml_file.stem.endswith("attmt"):
                continue
            found_value = search_xml_for_key(xml_file, search_key)
            existing_value = archive_metadata.get(search_key, None)
            if existing_value is None:
                archive_metadata[search_key] = found_value
            elif isinstance(existing_value, str):
                archive_metadata[search_key] = [existing_value, found_value]
            else:
                try:
                    archive_metadata[search_key].append(found_value)
                except Exception as e:
                    logger.error(e)
                    raise e

print(archive_metadata)
