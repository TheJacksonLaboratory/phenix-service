#!/usr/bin/env python3
import os
import json
import shutil
import sqlite3
import argparse
from pathlib import Path

import xmltodict


def construct_argparser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", title="subcommands")
    subparsers.required = True
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        "archive_root", type=Path,
        help="Path to 'Harmony-Archive' directory"
    )
    parent.add_argument(
        "-o", "--output-location", type=Path,
        help="Output directory to store images and/or metadata"
    )
    parent.add_argument(
        "--harmony-id", default=None,
        help="Optional harmony id string to select only one set of measurements"
    )

    convert = subparsers.add_parser(
        "convert",
        parents=[parent],
        help=(
            "Convert Harmony-Archive format to human-readable tiffs "
            "and metadata"
        )
    )
    convert.set_defaults(func=HarmonyArchive.convert_to_human_readable)

    metadata = subparsers.add_parser(
        "metadata",
        parents=[parent],
        help="Only generate metadata from a Harmony-Archive"
    )
    metadata.set_defaults(func=HarmonyArchive.generate_metadata_json)

    lister = subparsers.add_parser(
        "list",
        parents=[parent],
        help="List what exists in a Harmony Archive"
    )
    lister.set_defaults(output_location=None)
    lister.set_defaults(func=HarmonyArchive.list_experiments)

    return parser


class HarmonyArchive(object):
    human_readable_format = (
        "r{Row:02g}c{Col:02g}f{Field:02g}p{Plane:02g}-"
        "ch{Channel}sk{SlowKin}fk{FastKin}fl{Flim}.tiff"
    )

    xml_ns = {"http://www.perkinelmer.com/PEHH/HarmonyV5": None}

    def __init__(self, location):
        self.location = location
        self.exists = self.location.exists()
        if not self.exists:
            raise Exception(f"Harmony Archive directory {self.location} doesn't exist")
        self._validate_archive()

    def _validate_archive(self):
        # expects /IMAGES/IMAGES.sqlite
        # expects /XML/MEASUREMENT/*/*.xml
        self.image_db_locations = self.location.glob("IMAGES/*/IMAGES.sqlite")
        #assert self.image_db_location.exists(), \
        #    f"Image DB {self.image_db_location} not found!"
        self.measurement_xml_locations = self.location.glob("XML/MEASUREMENT/*.xml")

    def parse_xmls(self):
        xmls = {}
        for xml_loc in self.measurement_xml_locations:
            with open(xml_loc, "rb") as fin:
               x = xmltodict.parse(fin, process_namespaces=False, encoding="utf-8")#namespaces=self.xml_ns)
            root = x.get("Measurement", None)
            if root is None: continue
            mid = root.get("MeasurementID", None)
            if mid is not None:
                xmls[mid] = root
        self.metadata = xmls
        self.id_mapping = dict((k, v["PlateName"]) for k, v in xmls.items())

    def load_image_database(self, measurement=None):
        image_data = {}
        for database_location in self.image_db_locations:
            measurement_key = database_location.parent.name
            if (measurement is not None) and (measurement_key != measurement): continue

            with sqlite3.connect(str(database_location)) as image_db:
                image_db.row_factory = sqlite3.Row
                select_all = "SELECT * FROM Image"

                image_data[measurement_key] = [
                    dict(row) for row in image_db.execute(select_all)
                ]

        self.image_data = image_data

        self.parse_xmls()

    def convert_to_human_readable(self, output_location):
        if output_location is None:
            raise ValueError("Must specify output-location")
        if "image_data" not in self.__dict__:
            self.load_image_database()

        for key, records in self.image_data.items():
            for record in records:
                human_readable_image_name = self.human_readable_format.format(**record)
                record["human_readable"] = human_readable_image_name

                src_path = self.location / "IMAGES" / key / record["Url"]
                dest_name = self.id_mapping.get(key, key)
                dest_path = output_location / dest_name / human_readable_image_name

                if not dest_path.parent.exists():
                    os.makedirs(dest_path.parent)
                shutil.copyfile(src_path, dest_path)


    def generate_metadata_json(self, output_location):
        for xml_location in self.measurement_xml_locations:
            tree = ET.parse(xml_location)
            root = tree.getroot()
            #with open(xml_location, "r") as xml_in:
            for child in root.iter():
                #if child.tag.endswith("Measurement"):

                print(child.tag)
                break

    def list_experiments(self, *args, **kwargs):
        print(f"Measurements located in Archive: {self.location}")
        fmt = "{:<40}{:<30}{:>10}"
        print(fmt.format("Measurement ID", "Plate Name", "# images"))
        for key, records in self.image_data.items():
            print(fmt.format(key, self.id_mapping[key], len(records)))


def main():
    argparser = construct_argparser()
    args = argparser.parse_args()

    archive = HarmonyArchive(args.archive_root)
    archive.load_image_database(args.harmony_id)

    args.func(archive, args.output_location)


if __name__ == "__main__":
    main()
