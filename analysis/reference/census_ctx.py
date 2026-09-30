"""Open the CELLxGENE Census directly from its public S3 bucket (the mirror-discovery host is not reachable here)."""
import os
from urllib.parse import urlparse

import tiledbsoma as soma

CENSUS_URI = "s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/"


def open_census():
    cfg = {"vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"}
    proxy = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
    if proxy:
        p = urlparse(proxy)
        cfg.update({"vfs.s3.proxy_host": p.hostname, "vfs.s3.proxy_port": str(p.port), "vfs.s3.proxy_scheme": "http"})
        if os.path.exists("/root/.ccr/ca-bundle.crt"):
            cfg["vfs.s3.ca_file"] = "/root/.ccr/ca-bundle.crt"
    return soma.open(CENSUS_URI, context=soma.SOMATileDBContext(tiledb_config=cfg))
