#!/usr/bin/env python3
import os
import json
import shutil
import sqlite3
import argparse
from pathlib import Path

import xml.etree.ElementTree as ET


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
        "output_location", type=Path,
        help="Output directory to store images and/or metadata"
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
    return parser


class HarmonyArchive(object):
    human_readable_format = (
        "r{Row:02g}c{Col:02g}f{Field:02g}p{Plane:02g}-"
        "ch{Channel}sk{SlowKin}fk{FastKin}fl{Flim}.tiff"
    )

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
        self.measurement_xml_locations = self.location.glob("XML/*/*.xml")

    def load_image_database(self):
        image_data = {}
        for database_location in self.image_db_locations:
            measurement_key = database_location.parent.name

            with sqlite3.connect(str(database_location)) as image_db:
                image_db.row_factory = sqlite3.Row
                select_all = "SELECT * FROM Image"

                image_data[measurement_key] = [
                    dict(row) for row in image_db.execute(select_all)
                ]

        self.image_data = image_data

    def convert_to_human_readable(self, output_location):
        if "image_data" not in self.__dict__:
            self.load_image_database()

        for key, records in self.image_data.items():
            for record in records:
                human_readable_image_name = human_readable_format.format(**record)
                record["human_readable"] = human_readable_image_name

                src_path = self.location / "IMAGES" / key / record["Url"]
                dest_path = output_location / key / human_readable_image_name

                if not dest_path.parent.exists():
                    os.makedirs(dest_path.parent)
                shutil.copyfile(src_path, dest_path)
                break


    def generate_metadata_json(self, output_location):
        for xml_location in self.measurement_xml_locations:
            tree = ET.parse(xml_location)
            root = tree.getroot()
            #with open(xml_location, "r") as xml_in:
            for child in root.iter():
                #if child.tag.endswith("Measurement"):

                print(child.tag)
                break


def main():
    argparser = construct_argparser()
    args = argparser.parse_args()

    archive = HarmonyArchive(args.archive_root)
    archive.load_image_database()

    args.func(archive, args.output_location)


if __name__ == "__main__":
    main()
