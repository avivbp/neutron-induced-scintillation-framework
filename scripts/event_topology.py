"""Channel-aware neutron event topology helpers."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


CHANNELS = (
    "elastic",
    "inelastic",
    "capture",
    "fission",
    "other_hadronic",
    "transport_exit",
    "unclassified",
)
TOPOLOGY_CATEGORIES = (
    "no_recorded_interaction",
    "transport_exit_only",
    "single_elastic",
    "multiple_elastic",
    "inelastic",
    "elastic_plus_inelastic",
    "capture",
    "fission",
    "other_hadronic",
    "unclassified",
)


def _event_ids(values: pd.Series | list[int] | np.ndarray) -> pd.Index:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    return pd.Index(numeric.astype(np.int64).drop_duplicates(), name="event_id")


def classify_event_topologies(
    interactions: pd.DataFrame,
    event_ids: pd.Series | list[int] | np.ndarray | None = None,
) -> pd.DataFrame:
    """Aggregate interaction rows into one channel-aware row per event."""
    frame = interactions.copy()
    frame.columns = frame.columns.astype(str).str.strip()
    missing = sorted({"event_id", "channel"}.difference(frame.columns))
    if missing:
        raise ValueError(
            "interaction truth is missing columns: " + ", ".join(missing)
        )

    frame["event_id"] = pd.to_numeric(frame["event_id"], errors="coerce")
    frame = frame.loc[frame["event_id"].notna()].copy()
    frame["event_id"] = frame["event_id"].astype(np.int64)
    channel = frame["channel"].astype(str).str.strip().str.lower()
    channel = channel.where(channel.isin(CHANNELS), "unclassified")
    frame["channel"] = channel

    if event_ids is None:
        index = pd.Index(
            frame["event_id"].drop_duplicates().sort_values(), name="event_id"
        )
    else:
        index = _event_ids(event_ids)

    counts = pd.crosstab(frame["event_id"], frame["channel"])
    counts = counts.reindex(index=index, columns=CHANNELS, fill_value=0)
    result = pd.DataFrame(index=index)
    for channel_name in CHANNELS:
        result[f"{channel_name}_count"] = counts[channel_name].astype(np.int64)

    roles = (
        frame["volume_roles"].fillna("").astype(str)
        if "volume_roles" in frame.columns
        else pd.Series("", index=frame.index)
    )
    fiducial_elastic = frame.loc[
        frame["channel"].eq("elastic")
        & roles.str.contains(r"(?:^|\|)fiducial_lar(?:\||$)", regex=True),
        "event_id",
    ].value_counts()
    result["fiducial_elastic_count"] = (
        fiducial_elastic.reindex(index, fill_value=0).astype(np.int64)
    )

    if "parent_id" in frame.columns:
        parent_id = pd.to_numeric(frame["parent_id"], errors="coerce")
        primary = frame.loc[parent_id.eq(0), "event_id"].value_counts()
        secondary = frame.loc[parent_id.gt(0), "event_id"].value_counts()
        result["primary_neutron_interaction_count"] = primary.reindex(
            index, fill_value=0
        ).astype(np.int64)
        result["secondary_neutron_interaction_count"] = secondary.reindex(
            index, fill_value=0
        ).astype(np.int64)
    else:
        result["primary_neutron_interaction_count"] = 0
        result["secondary_neutron_interaction_count"] = 0

    elastic = result["elastic_count"]
    inelastic = result["inelastic_count"]
    recorded = result[[f"{value}_count" for value in CHANNELS]].sum(axis=1)
    result["event_type"] = np.select(
        (
            result["fission_count"].gt(0),
            result["capture_count"].gt(0),
            result["other_hadronic_count"].gt(0),
            result["unclassified_count"].gt(0),
            inelastic.gt(0) & elastic.gt(0),
            inelastic.gt(0),
            elastic.gt(1),
            elastic.eq(1),
            result["transport_exit_count"].gt(0),
            recorded.eq(0),
        ),
        (
            "fission",
            "capture",
            "other_hadronic",
            "unclassified",
            "elastic_plus_inelastic",
            "inelastic",
            "multiple_elastic",
            "single_elastic",
            "transport_exit_only",
            "no_recorded_interaction",
        ),
        default="unclassified",
    )
    return result.reset_index()


def _safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", value)


def interaction_truth_path_for_pe(pe_path: Path) -> Path | None:
    """Resolve the interaction-truth file matching one numPE scan file."""
    match = re.search(
        r"numPE_(\d+)_topPMTs_(\d+)_botPMTs_(\d+)_SiPMRows",
        pe_path.name,
    )
    metadata_path = pe_path.parent / "scan_metadata.json"
    if match and metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            configuration = tuple(int(value) for value in match.groups())
            matches = [
                row
                for row in metadata
                if (
                    int(row["top_pmts"]),
                    int(row["bottom_pmts"]),
                    int(row["sipm_rows"]),
                )
                == configuration
            ]
            if len(matches) == 1:
                candidate = pe_path.parent / (
                    "neutron_interactions_"
                    + _safe_label(str(matches[0]["config_id"]))
                    + ".csv"
                )
                if candidate.exists():
                    return candidate
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
            pass

    candidates = sorted(pe_path.parent.glob("neutron_interactions_*.csv"))
    pe_candidates = list(
        pe_path.parent.glob("numPE_*_topPMTs_*_botPMTs_*_SiPMRows.csv")
    )
    if len(candidates) == 1 and len(pe_candidates) == 1:
        return candidates[0]
    return None


def true_single_elastic_mask_from_truth(
    events: pd.DataFrame, interaction_path: Path
) -> pd.Series:
    """Select events with one fiducial elastic and no other hadronic channel."""
    if "eventID" not in events.columns:
        raise ValueError("interaction-truth selection requires eventID")
    interactions = pd.read_csv(interaction_path)
    event_ids = pd.to_numeric(events["eventID"], errors="coerce")
    topology = classify_event_topologies(interactions, event_ids)
    selected_ids = topology.loc[
        topology["event_type"].eq("single_elastic")
        & topology["fiducial_elastic_count"].eq(1),
        "event_id",
    ]
    return event_ids.isin(selected_ids)
