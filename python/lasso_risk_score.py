# ==============================================================================
# lasso_risk_score.py
#
# LASSO dementia-risk score (Figure S22; predicted_cog_vulnerability /
# high_baseline_vulnerability). Faithful port of the project's original
# lasso_predict_dementiarisk.ipynb: same library (sklearn ElasticNetCV,
# l1_ratio=1), same pipeline (median-impute -> standardize -> CV), same
# random_state. R's glmnet gives a different, much less sparse solution on
# this dataset (cv.glmnet's random K-fold assignment vs. sklearn's default
# unshuffled KFold, on a weak-signal target where the CV error curve is
# flat enough that this choice changes which features survive) -- so this
# step runs in Python rather than being re-implemented in R.
#
# Called from R/00c_lasso_risk_score.R, which exports the pre-vaccination
# (+ earliest-visit-fallback) panel to data/_lasso_long.csv and the
# regression target to data/_lasso_target.csv, then reads back this
# script's data/_lasso_coefficients.csv and data/_lasso_predictions.csv.
# ==============================================================================
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.impute import SimpleImputer

import random
random.seed(100)
np.random.seed(100)

d = pd.read_csv("data/_lasso_long.csv")
d["visit_date"] = pd.to_datetime(d["visit_date"])
target = pd.read_csv("data/_lasso_target.csv")

core_features = [
    "b_visit_age", "female", "white", "education",
    "b_casi_total", "apoe_any4",
    "firstdeg_dem", "firstdeg_ad", "firstdeg_memprob",
    "income", "basicneeds_deprived",
    "lb_depre", "cesd20_score",
    "hypertension", "diabetes", "stroke", "mi", "chf", "sleepapnea",
    "avesys", "avedia", "bmi",
    "smoke", "smoke_intensity", "alcohol", "exercise_reg",
    "adl_sum", "iadl_sum", "gait_time", "grip_strength",
]
immune_features = ["rheumarth", "lupus", "othautoimm", "organbmt", "chemo", "othimmsup"]
hz_features = ["hz_beforevax"]
ccw_features = [c for c in d.columns if c.startswith("ccw_")]
ip_features = [c for c in d.columns if c.startswith("ip_")]

feature_list = [v for v in (core_features + immune_features + hz_features + ccw_features + ip_features)
                 if v in d.columns]
missing_rate = d[feature_list].isna().mean()
feature_list = [v for v in feature_list if missing_rate[v] < 1]
print("N features used:", len(feature_list))

first_vars = ["b_visit_age", "female", "white", "education", "b_casi_total", "apoe_any4",
              "firstdeg_dem", "firstdeg_ad", "firstdeg_memprob"]
mean_vars = ["income", "basicneeds_deprived", "cesd20_score", "avesys", "avedia", "bmi",
             "adl_sum", "iadl_sum", "gait_time", "grip_strength", "smoke_intensity"]
last_vars = (["lb_depre", "smoke", "alcohol", "exercise_reg", "hz_beforevax",
              "lupus", "othautoimm", "organbmt", "chemo", "othimmsup"]
             + ccw_features + ip_features)

d = d.sort_values("visit_date")
base_first = d.groupby("subject")[[v for v in first_vars if v in d.columns]].first()
base_mean = d.groupby("subject")[[v for v in mean_vars if v in d.columns]].mean(numeric_only=True)
base_last = d.groupby("subject")[[v for v in last_vars if v in d.columns]].last()

baseline_all = pd.concat([base_first, base_mean, base_last], axis=1).reset_index()
baseline_all = baseline_all.merge(target, on="subject", how="left")

X_all = baseline_all[[v for v in feature_list if v in baseline_all.columns]].copy()
y_all = baseline_all["prior_cog_vulnerability"].copy()

obs = y_all.notna()
X_obs, y_obs = X_all.loc[obs].copy(), y_all.loc[obs].copy()
X_tr, X_te, y_tr, y_te = train_test_split(X_obs, y_obs, test_size=0.2, random_state=100)

fill_vals = X_tr.median(numeric_only=True)
X_tr = X_tr.fillna(fill_vals)
X_te = X_te.fillna(fill_vals)
X_all = X_all.fillna(fill_vals)

enet = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
    ("enet", ElasticNetCV(l1_ratio=1, cv=5, random_state=100, n_jobs=1)),
])
enet.fit(X_tr, y_tr)

yhat = enet.predict(X_te)
print("TEST RMSE:", round(float(np.sqrt(mean_squared_error(y_te, yhat))), 4))
print("TEST R2:", round(float(r2_score(y_te, yhat)), 4))

coef = enet.named_steps["enet"].coef_
importance = pd.Series(coef, index=X_tr.columns)
importance = importance[importance != 0].sort_values(key=np.abs, ascending=False)
print("\nNon-zero LASSO coefficients (standardized X), largest |coef| first:")
print(importance)
importance.rename_axis("term").reset_index(name="estimate").to_csv(
    "data/_lasso_coefficients.csv", index=False)

# ------------------------------------------------------------------------
# Figure S22: same matplotlib styling and predictor labels as the original
# notebook (rename_map / auto_label), not an R/ggplot re-creation.
# ------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "serif"],
    "font.size": 11,
    "axes.titlesize": 15,
    "xtick.labelsize": 14,
    "axes.labelsize": 15,
    "ytick.labelsize": 15,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

rename_map = {
    "b_visit_age": "Age at enrollment", "b_casi_total": "CASI score at enrollment",
    "white": "White (enrollment)", "education": "Education (enrollment)",
    "adl_sum": "ADL limitations (avg)", "iadl_sum": "IADL limitations (avg)",
    "gait_time": "Gait time (avg)", "grip_strength": "Grip strength (avg)",
    "cesd20_score": "Depression score (CESD, avg)", "bmi": "BMI (avg)",
    "avesys": "Systolic BP (avg)", "avedia": "Diastolic BP (avg)",
    "smoke_intensity": "Smoking intensity (avg)", "exercise_reg": "Regular exercise (avg)",
    "smoke": "Smoking (last pre-vax)", "alcohol": "Alcohol use (last pre-vax)",
    "hypertension": "Hypertension (self-report, last pre-vax)",
    "diabetes": "Diabetes (self-report, last pre-vax)",
    "stroke": "Stroke (self-report, last pre-vax)",
    "mi": "Myocardial infarction (self-report, last pre-vax)",
    "chf": "Heart failure (self-report, last pre-vax)",
    "sleepapnea": "Sleep apnea (self-report, last pre-vax)",
    "rheumarth": "Rheumatoid arthritis (self-report, last pre-vax)",
    "lb_depre": "Depression (diagnosis, last pre-vax)",
    "hz_beforevax": "Prior herpes zoster (last pre-vax)",
    "ccw_hthy": "Hypothyroidism (diagnosis)", "ccw_depr": "Depression (diagnosis)",
    "ccw_canc": "Cancer (diagnosis)", "ccw_myoc": "Myocardial infarction (diagnosis)",
    "ccw_raoa": "Rheumatoid arthritis / Osteoarthritis (diagnosis)",
    "ccw_stro": "Stroke (diagnosis)", "ccw_htn": "Hypertension (diagnosis)",
    "ip_septi": "Sepsis (inpatient diagnosis)",
}

def auto_label(v):
    if v in rename_map:
        return rename_map[v]
    if v.startswith("ccw_"):
        return v.replace("ccw_", "").replace("_", " ").title() + " (diagnosis)"
    if v.startswith("ip_"):
        return v.replace("ip_", "").replace("_", " ").title() + " (inpatient diagnosis)"
    return v

plot_importance = importance.copy()
plot_importance.index = [auto_label(v) for v in plot_importance.index]

fig, ax = plt.subplots(figsize=(14, 8))
colors = ["#d62728" if v < 0 else "#1f77b4" for v in plot_importance.values]
ax.barh(plot_importance.index, plot_importance.values, color=colors)
ax.invert_yaxis()
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Coefficient value")
ax.tick_params(axis="both")
plt.tight_layout()
plt.savefig("output/figures/FigureS22_lasso_coefficients.png", dpi=600, bbox_inches="tight")
print("Saved output/figures/FigureS22_lasso_coefficients.png")

baseline_all["predicted_cog_vulnerability"] = enet.predict(X_all)
baseline_all[["subject", "predicted_cog_vulnerability"]].to_csv(
    "data/_lasso_predictions.csv", index=False)
