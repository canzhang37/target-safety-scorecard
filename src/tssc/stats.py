"""AUROC with bootstrap, grouped cross-validation, permutation tests."""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def boot_auc(y, scores, reps=1000, seed=1):
    """Point AUROC and bootstrap draws for several scores, same resamples."""
    y = np.asarray(y)
    rng = np.random.default_rng(seed)
    point = {k: roc_auc_score(y, s) for k, s in scores.items()}
    draws = {k: [] for k in scores}
    n = len(y)
    for _ in range(reps):
        i = rng.integers(0, n, n)
        if y[i].min() == y[i].max():
            continue
        for k, s in scores.items():
            draws[k].append(roc_auc_score(y[i], np.asarray(s)[i]))
    return point, {k: np.array(v) for k, v in draws.items()}


def ci(d):
    return np.percentile(d, [2.5, 97.5])


def oof_logistic(X, y, groups, folds=5):
    """Out-of-fold probabilities, folds grouped by target class."""
    X, y = np.asarray(X, float), np.asarray(y)
    p = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=folds).split(X, y, groups):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        m.fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p


def impute_with_flags(df, cols):
    """Median-impute and add a missing flag per column with gaps."""
    out = df[cols].copy()
    for c in cols:
        miss = out[c].isna()
        if miss.any():
            out[c + "_missing"] = miss.astype(float)
            out[c] = out[c].fillna(out[c].median())
    return out


def concordant_pairs(labels, groups):
    """Same-group pairs with equal labels, and all same-group pairs."""
    labels, groups = np.asarray(labels), np.asarray(groups)
    _, gi = np.unique(groups, return_inverse=True)
    k = np.bincount(gi)
    m = np.bincount(gi, weights=labels)
    same = (m * (m - 1) + (k - m) * (k - m - 1)) / 2
    return same.sum(), (k * (k - 1) / 2).sum()


def concordance_test(labels, groups, perms=10000, seed=1):
    labels = np.asarray(labels, float)
    obs, tot = concordant_pairs(labels, groups)
    rng = np.random.default_rng(seed)
    null = np.array([concordant_pairs(rng.permutation(labels), groups)[0] for _ in range(perms)])
    p = (1 + (null >= obs).sum()) / (perms + 1)
    return obs / tot, null.mean() / tot, p, int(tot)


def stratified_diff_test(z, member, strata, perms=10000, seed=1):
    """Mean of within-stratum differences (members minus others),
    weighted by member count; labels shuffled within strata."""
    z, member, strata = np.asarray(z, float), np.asarray(member, bool), np.asarray(strata)
    levels = [s for s in np.unique(strata)
              if member[strata == s].any() and (~member[strata == s]).any()]

    def stat(mem):
        num = den = 0.0
        for s in levels:
            i = strata == s
            w = mem[i].sum()
            num += w * (z[i][mem[i]].mean() - z[i][~mem[i]].mean())
            den += w
        return num / den

    obs = stat(member)
    rng = np.random.default_rng(seed)
    idx = {s: np.flatnonzero(strata == s) for s in levels}
    null = np.empty(perms)
    for b in range(perms):
        mem = member.copy()
        for s, i in idx.items():
            mem[i] = rng.permutation(member[i])
        null[b] = stat(mem)
    return obs, (1 + (null >= obs).sum()) / (perms + 1), levels
