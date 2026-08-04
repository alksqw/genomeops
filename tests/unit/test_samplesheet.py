from __future__ import annotations

import csv
from pathlib import Path

import pytest

from genomeops.samplesheet import (
    REQUIRED_COLUMNS,
    SamplesheetValidationError,
    validate_samplesheet,
    write_normalized_samplesheet,
)


def write_test_samplesheet(
    path: Path,
    rows: list[dict[str, str]],
    *,
    fieldnames: list[str] | tuple[str, ...] = REQUIRED_COLUMNS,
) -> None:
    with path.open(
        mode="w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def valid_row() -> dict[str, str]:
    return {
        "sample": "HG002",
        "lane": "L001",
        "fastq_1": "data/HG002_L001_R1.fastq.gz",
        "fastq_2": "data/HG002_L001_R2.fastq.gz",
        "library": "LIB001",
        "platform": "illumina",
        "center": "ITMO",
    }


def test_valid_samplesheet_is_normalized(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(input_path, [valid_row()])

    rows = validate_samplesheet(input_path)

    assert len(rows) == 1
    assert rows[0]["sample"] == "HG002"
    assert rows[0]["platform"] == "ILLUMINA"
    assert "ID:HG002.L001" in rows[0]["read_group"]
    assert "SM:HG002" in rows[0]["read_group"]


def test_multiple_lanes_for_one_sample_are_allowed(
    tmp_path: Path,
) -> None:
    first_row = valid_row()
    second_row = valid_row()

    second_row["lane"] = "L002"
    second_row["fastq_1"] = "data/HG002_L002_R1.fastq.gz"
    second_row["fastq_2"] = "data/HG002_L002_R2.fastq.gz"

    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(
        input_path,
        [first_row, second_row],
    )

    rows = validate_samplesheet(input_path)

    assert len(rows) == 2
    assert rows[0]["lane"] == "L001"
    assert rows[1]["lane"] == "L002"


def test_missing_required_column_is_rejected(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "samplesheet.csv"

    fieldnames = [
        column
        for column in REQUIRED_COLUMNS
        if column != "center"
    ]

    row = valid_row()
    row.pop("center")

    write_test_samplesheet(
        input_path,
        [row],
        fieldnames=fieldnames,
    )

    with pytest.raises(
        SamplesheetValidationError,
        match="Missing required column: center",
    ):
        validate_samplesheet(input_path)


def test_duplicate_sample_lane_is_rejected(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "samplesheet.csv"

    write_test_samplesheet(
        input_path,
        [valid_row(), valid_row()],
    )

    with pytest.raises(
        SamplesheetValidationError,
        match="duplicate sample/lane combination",
    ):
        validate_samplesheet(input_path)


def test_same_fastq_is_rejected(
    tmp_path: Path,
) -> None:
    row = valid_row()
    row["fastq_2"] = row["fastq_1"]

    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(input_path, [row])

    with pytest.raises(
        SamplesheetValidationError,
        match="fastq_1 and fastq_2",
    ):
        validate_samplesheet(input_path)


def test_invalid_sample_identifier_is_rejected(
    tmp_path: Path,
) -> None:
    row = valid_row()
    row["sample"] = "HG002 bad sample"

    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(input_path, [row])

    with pytest.raises(
        SamplesheetValidationError,
        match="invalid sample identifier",
    ):
        validate_samplesheet(input_path)


def test_unsupported_platform_is_rejected(
    tmp_path: Path,
) -> None:
    row = valid_row()
    row["platform"] = "UNKNOWN"

    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(input_path, [row])

    with pytest.raises(
        SamplesheetValidationError,
        match="unsupported platform",
    ):
        validate_samplesheet(input_path)


def test_missing_local_files_are_rejected_when_requested(
    tmp_path: Path,
) -> None:
    row = valid_row()
    row["fastq_1"] = str(tmp_path / "missing_R1.fastq.gz")
    row["fastq_2"] = str(tmp_path / "missing_R2.fastq.gz")

    input_path = tmp_path / "samplesheet.csv"
    write_test_samplesheet(input_path, [row])

    with pytest.raises(
        SamplesheetValidationError,
        match="local file does not exist",
    ):
        validate_samplesheet(
            input_path,
            check_files=True,
        )


def test_normalized_samplesheet_is_written(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "samplesheet.csv"
    output_path = tmp_path / "normalized.csv"

    write_test_samplesheet(input_path, [valid_row()])

    rows = validate_samplesheet(input_path)
    write_normalized_samplesheet(rows, output_path)

    assert output_path.is_file()

    output_text = output_path.read_text(encoding="utf-8")

    assert "read_group" in output_text
    assert "HG002" in output_text
    assert "ILLUMINA" in output_text
