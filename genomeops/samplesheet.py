from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from urllib.parse import urlparse


REQUIRED_COLUMNS = (
    "sample",
    "lane",
    "fastq_1",
    "fastq_2",
    "library",
    "platform",
    "center",
)

OUTPUT_COLUMNS = REQUIRED_COLUMNS + ("read_group",)

ALLOWED_PLATFORMS = {"ILLUMINA"}

FASTQ_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)

SAFE_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

REMOTE_SCHEMES = {
    "http",
    "https",
    "ftp",
    "s3",
    "gs",
    "az",
}


class SamplesheetValidationError(ValueError):
    """Raised when a samplesheet does not satisfy the required format."""

    def __init__(self, errors: list[str]) -> None:
        self.errors = errors
        super().__init__("\n".join(errors))


def _is_remote_location(value: str) -> bool:
    """Return True when a path points to supported remote storage."""

    scheme = urlparse(value).scheme.lower()
    return scheme in REMOTE_SCHEMES


def _has_fastq_suffix(value: str) -> bool:
    """Return True when the location has a supported FASTQ extension."""

    parsed_path = urlparse(value).path.lower()
    return parsed_path.endswith(FASTQ_SUFFIXES)


def _create_read_group(
    *,
    sample: str,
    lane: str,
    library: str,
    platform: str,
    center: str,
) -> str:
    """Create a BWA-compatible read-group string."""

    read_group_id = f"{sample}.{lane}"

    return (
        f"@RG\\t"
        f"ID:{read_group_id}\\t"
        f"SM:{sample}\\t"
        f"LB:{library}\\t"
        f"PL:{platform}\\t"
        f"PU:{read_group_id}\\t"
        f"CN:{center}"
    )


def validate_samplesheet(
    input_path: Path,
    *,
    check_files: bool = False,
) -> list[dict[str, str]]:
    """
    Validate a GenomeOps samplesheet and return normalized rows.

    Structural validation is always performed. Local FASTQ existence is
    checked only when check_files=True.
    """

    errors: list[str] = []

    if not input_path.is_file():
        raise SamplesheetValidationError(
            [f"Input samplesheet does not exist: {input_path}"]
        )

    normalized_rows: list[dict[str, str]] = []
    seen_sample_lanes: set[tuple[str, str]] = set()
    nonempty_row_count = 0

    with input_path.open(
        mode="r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise SamplesheetValidationError(
                ["The samplesheet is empty or has no header."]
            )

        normalized_headers = [
            header.strip() if header is not None else ""
            for header in reader.fieldnames
        ]
        reader.fieldnames = normalized_headers

        duplicate_headers = sorted(
            {
                header
                for header in normalized_headers
                if normalized_headers.count(header) > 1
            }
        )

        for header in duplicate_headers:
            errors.append(f"Duplicate column in header: {header}")

        missing_columns = [
            column
            for column in REQUIRED_COLUMNS
            if column not in normalized_headers
        ]

        for column in missing_columns:
            errors.append(f"Missing required column: {column}")

        if errors:
            raise SamplesheetValidationError(errors)

        for row_number, raw_row in enumerate(reader, start=2):
            extra_values = raw_row.get(None)

            if extra_values and any(
                str(value).strip() for value in extra_values
            ):
                errors.append(
                    f"Row {row_number}: contains more values than the header."
                )

            row = {
                key: (value or "").strip()
                for key, value in raw_row.items()
                if key is not None
            }

            if not any(row.get(column, "") for column in REQUIRED_COLUMNS):
                continue

            nonempty_row_count += 1
            error_count_before_row = len(errors)

            empty_columns = [
                column
                for column in REQUIRED_COLUMNS
                if not row.get(column, "")
            ]

            for column in empty_columns:
                errors.append(
                    f"Row {row_number}: required value is empty: {column}"
                )

            sample = row.get("sample", "")
            lane = row.get("lane", "")
            fastq_1 = row.get("fastq_1", "")
            fastq_2 = row.get("fastq_2", "")
            library = row.get("library", "")
            platform = row.get("platform", "").upper()
            center = row.get("center", "")

            if sample and not SAFE_IDENTIFIER_PATTERN.fullmatch(sample):
                errors.append(
                    f"Row {row_number}: invalid sample identifier: {sample}"
                )

            if lane and not SAFE_IDENTIFIER_PATTERN.fullmatch(lane):
                errors.append(
                    f"Row {row_number}: invalid lane identifier: {lane}"
                )

            if platform and platform not in ALLOWED_PLATFORMS:
                allowed = ", ".join(sorted(ALLOWED_PLATFORMS))
                errors.append(
                    f"Row {row_number}: unsupported platform "
                    f"'{platform}'. Allowed: {allowed}"
                )

            for field_name, location in (
                ("fastq_1", fastq_1),
                ("fastq_2", fastq_2),
            ):
                if not location:
                    continue

                if not _has_fastq_suffix(location):
                    errors.append(
                        f"Row {row_number}: {field_name} has an unsupported "
                        f"FASTQ extension: {location}"
                    )

                if (
                    check_files
                    and not _is_remote_location(location)
                    and not Path(location).is_file()
                ):
                    errors.append(
                        f"Row {row_number}: local file does not exist "
                        f"for {field_name}: {location}"
                    )

            if fastq_1 and fastq_2 and fastq_1 == fastq_2:
                errors.append(
                    f"Row {row_number}: fastq_1 and fastq_2 "
                    "refer to the same location."
                )

            if sample and lane:
                sample_lane_key = (sample, lane)

                if sample_lane_key in seen_sample_lanes:
                    errors.append(
                        f"Row {row_number}: duplicate sample/lane "
                        f"combination: {sample}/{lane}"
                    )
                else:
                    seen_sample_lanes.add(sample_lane_key)

            if len(errors) == error_count_before_row:
                normalized_rows.append(
                    {
                        "sample": sample,
                        "lane": lane,
                        "fastq_1": fastq_1,
                        "fastq_2": fastq_2,
                        "library": library,
                        "platform": platform,
                        "center": center,
                        "read_group": _create_read_group(
                            sample=sample,
                            lane=lane,
                            library=library,
                            platform=platform,
                            center=center,
                        ),
                    }
                )

    if nonempty_row_count == 0:
        errors.append("The samplesheet contains no data rows.")

    if errors:
        raise SamplesheetValidationError(errors)

    return normalized_rows


def write_normalized_samplesheet(
    rows: list[dict[str, str]],
    output_path: Path,
) -> None:
    """Write validated rows to a normalized CSV file."""

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        mode="w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Validate and normalize a GenomeOps FASTQ samplesheet."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to the input CSV samplesheet.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path for the normalized output CSV.",
    )

    parser.add_argument(
        "--check-files",
        action="store_true",
        help=(
            "Require local FASTQ paths to exist. "
            "Remote URLs are not downloaded."
        ),
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the samplesheet validator CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        rows = validate_samplesheet(
            args.input,
            check_files=args.check_files,
        )
    except SamplesheetValidationError as error:
        print(
            "[ERROR] Samplesheet validation failed:",
            file=sys.stderr,
        )

        for message in error.errors:
            print(f"  - {message}", file=sys.stderr)

        return 1

    write_normalized_samplesheet(rows, args.output)

    print(f"[OK] Validated {len(rows)} samplesheet row(s).")
    print(f"[OK] Normalized samplesheet: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
