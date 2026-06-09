from __future__ import annotations

# Load environment variables from .env file
from dotenv import load_dotenv
from pathlib import Path
env_path = Path(__file__).parent.parent / ".env"
load_dotenv(dotenv_path=env_path)


from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from .pcr import router as pcr_router
from .gibson_plan import router as gibson_plan_router
from .gibson_primers import router as gibson_primers_router
from .batch_pcr import router as batch_pcr_router
from .plannotate_router import router as plannotate_router
from .inv_gib import router as inv_gib_router
from .chat import router as chat_router
from .cloning.router import router as cloning_router
from .sbol_router import router as sbol_router
try:
    from .plasmid_lm.lm_router import router as plasmid_lm_router
except Exception:
    plasmid_lm_router = None

try:
    from .agent.agent_router import router as agent_router
except Exception:
    agent_router = None

try:
    from .biosecurity_router import router as biosecurity_router
except ImportError:
    biosecurity_router = None

try:
    from .deploy_router import router as deploy_router
except ImportError:
    deploy_router = None

app = FastAPI(title="Splicify API")

# CORS middleware for frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Sanity endpoint (kept) ---
class TextRequest(BaseModel):
    text: str

@app.post("/process")
def process(req: TextRequest):
    return {"length": len(req.text), "uppercase": req.text.upper()}

# --- Include feature routers ---
app.include_router(chat_router)  # Main chat orchestration endpoint
app.include_router(batch_pcr_router)
app.include_router(pcr_router)
app.include_router(gibson_plan_router)
app.include_router(gibson_primers_router)
app.include_router(plannotate_router)
app.include_router(inv_gib_router)
app.include_router(cloning_router)
app.include_router(sbol_router)
if plasmid_lm_router is not None:
    app.include_router(plasmid_lm_router)
if agent_router is not None:
    app.include_router(agent_router)
if biosecurity_router is not None:
    app.include_router(biosecurity_router)
if deploy_router is not None:
    app.include_router(deploy_router)

# Pre-load all annotation components at startup to avoid first-request delay
@app.on_event("startup")
async def startup_event():
    try:
        # Pre-load KnowledgeBase
        from . import plasmid_analyzer
        from .plasmid_analyzer import KnowledgeBase, build_feature_instances, detect_modules
        KnowledgeBase._load_kb()
        print("KnowledgeBase pre-loaded at startup")
        
        # Pre-load hierarchical annotator rules
        from .hierarchical_annotator import _load_rules
        _load_rules()
        print("Hierarchical annotator rules pre-loaded")
        
        # Warm up the module extractor imports
        from .Module_Library_gb import module_extractor
        print("Module extractor pre-loaded")
        
        # Warmup: run a minimal analysis to trigger all lazy initialization
        from .hierarchical_annotator import annotate_hierarchy_from_plannotate_v2
        warmup_result = annotate_hierarchy_from_plannotate_v2(
            sequence="ATGC",
            circular=True,
            plannotate_rows=[]
        )
        print("Warmup analysis completed")
    except Exception as e:
        print(f"Warning: Could not pre-load components: {e}")
