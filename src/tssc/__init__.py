"""Target safety scorecard on Open Targets data."""
from pathlib import Path

RELEASE = "26.09"
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / RELEASE
RES = ROOT / "results"

KEEP_DRUG_TYPES = {"Small molecule", "Antibody", "Protein", "Oligonucleotide", "Oligosaccharide"}
DROP_ACTIONS = {"EXOGENOUS GENE", "EXOGENOUS PROTEIN", "VACCINE ANTIGEN", "HYDROLYTIC ENZYME",
                "PROTEOLYTIC ENZYME", "CROSS-LINKING AGENT", "DISRUPTING AGENT"}
WARNING_TYPES = {"Black Box Warning", "Withdrawn"}

# risk = minus the Open Targets factor
PRE_DRUG = ["geneticConstraint", "mouseKOScore", "isCancerDriverGene",
            "paralogMaxIdentityPercentage", "tissueSpecificity", "tissueDistribution"]

CLINICAL_SOURCES = {"Bowes et al. (2012)", "Brennan et al. (2024)", "Urban et al. (2012)",
                    "Lynch et al. (2017)", "Force et al. (2011)", "ClinPGx"}
EXPERIMENTAL_SOURCES = {"ToxCast", "AOP-Wiki", "Lamore et al. (2017)"}

# toxicity class -> keywords; a keyword equal to the
# parent organ or inside the tissue name matches
ORGANS = {
    "hepatotoxicity": ["liver"],
    "cardiotoxicity": ["heart"],
    "nephrotoxicity": ["kidney"],
    "neurotoxicity": ["brain"],
    "psychiatric toxicity": ["brain"],
    "dermatological toxicity": ["skin"],
    "gastrointestinal toxicity": ["colon", "stomach", "esophagus", "small intestine"],
    "respiratory toxicity": ["lung"],
    "hematological toxicity": ["blood"],
    "musculoskeletal toxicity": ["muscle tissue"],
    "vascular toxicity": ["artery", "aorta"],
}
