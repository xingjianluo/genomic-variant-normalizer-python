import uvicorn
from fastapi import FastAPI

from variant_service.api import router as variant_router

app = FastAPI(title="Offline Genomic Variant Normalizer")
app.include_router(variant_router)


@app.get("/")
def read_root():
    return {"service": "variant-normalizer", "status": "ok"}


if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.1", port=0)
