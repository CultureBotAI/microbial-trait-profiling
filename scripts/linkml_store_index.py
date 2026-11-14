#!/usr/bin/env python3
"""Build a search index over the bsdbng studies collection in MongoDB.

Requires data loaded by linkml_store_load.py.

Indexes are persisted in MongoDB as `internal__index__studies__{kind}`,
so this is a one-time operation per kind. Re-running with --rebuild drops
and recreates the index.

Kinds:
    simple   -- trigram text index (no API key required, ~seconds)
    llm      -- OpenAI text-embedding-ada-002 semantic index
                (requires OPENAI_API_KEY, ~minutes for 1,936 studies)

Usage:
    uv run python scripts/linkml_store_index.py --kind simple
    uv run python scripts/linkml_store_index.py --kind llm
    uv run python scripts/linkml_store_index.py --kind llm --rebuild
"""

from __future__ import annotations

import argparse
import os
import time

from dotenv import load_dotenv
from linkml_store import Client
from linkml_store.index import get_indexer

DEFAULT_MONGO_URI = "mongodb://localhost:27017/bsdbng"
COLLECTION_NAME = "studies"


def main() -> None:
    load_dotenv()

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--kind",
        choices=["simple", "llm"],
        required=True,
        help="Index kind to build",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Drop and rebuild the index even if already present",
    )
    args = parser.parse_args()

    mongo_uri = os.environ.get("BSDBNG_MONGO_URI", DEFAULT_MONGO_URI)

    client = Client()
    db = client.attach_database(mongo_uri, alias="bsdbng")
    collection = db.get_collection(COLLECTION_NAME)

    from pymongo import MongoClient

    mongo: MongoClient[dict[str, object]] = MongoClient(mongo_uri)
    mdb = mongo.get_database()
    index_coll_name = f"internal__index__{COLLECTION_NAME}__{args.kind}"
    existing = index_coll_name in mdb.list_collection_names()

    if existing and args.rebuild:
        mdb.drop_collection(index_coll_name)
        mdb.drop_collection(index_coll_name + "__metadata")
        print(f"Dropped existing {args.kind} index")
        existing = False

    if existing:
        n = mdb.get_collection(index_coll_name).count_documents({})
        print(f"{args.kind} index already has {n} vectors; use --rebuild to recreate")
        return

    t0 = time.time()
    collection.attach_indexer(get_indexer(args.kind), args.kind)
    print(f"Built {args.kind} index in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
