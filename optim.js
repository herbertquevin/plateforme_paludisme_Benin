// Optimisation budgétaire (remplace scipy.optimize.minimize SLSQP du backend
// Python) par une descente de gradient projetée : le problème (minimiser un
// rayon spectral, lisse et de signe constant de gradient, sous 0<=u_k<=0.95 et
// somme(u)<=budget) s'y prête bien et évite d'embarquer un solveur NLP complet
// dans la page. Vérifié par comparaison directe aux résultats SLSQP de
// scipy (voir gen_reference.py) : écarts < 1e-2 sur R0 optimisé et
// allocation, bien en-deçà de la précision utile pour la décision.

function projectBoxSum(v, lo, hi, budget) {
  const n = v.length;
  const clip = (x) => Math.max(lo, Math.min(hi, x));
  const clipped = v.map(clip);
  const s = clipped.reduce((a, b) => a + b, 0);
  if (s <= budget + 1e-12) return clipped;
  // bisection sur tau : sum(clip(v - tau, lo, hi)) = budget, décroissant en tau
  let lo_t = -1.0, hi_t = 1.0;
  // élargit la borne si besoin
  while (v.map((x) => clip(x - lo_t)).reduce((a, b) => a + b, 0) < budget) lo_t -= 1.0;
  while (v.map((x) => clip(x - hi_t)).reduce((a, b) => a + b, 0) > budget) hi_t += 1.0;
  for (let it = 0; it < 100; it++) {
    const mid = (lo_t + hi_t) / 2;
    const sMid = v.map((x) => clip(x - mid)).reduce((a, b) => a + b, 0);
    if (sMid > budget) lo_t = mid; else hi_t = mid;
  }
  const tau = (lo_t + hi_t) / 2;
  return v.map((x) => clip(x - tau));
}

function optimiserAllocation(f, n, budget, opts = {}) {
  const lo = 0.0, hi = 0.95;
  const maxIter = opts.maxIter || 300;
  const h = opts.h || 1e-4;
  let u = projectBoxSum(new Array(n).fill(budget / n), lo, hi, budget);
  let fu = f(u);
  let step = opts.step0 || 0.5;

  for (let it = 0; it < maxIter; it++) {
    // gradient par différences finies centrées
    const grad = new Array(n);
    for (let k = 0; k < n; k++) {
      const up = u.slice(); up[k] = Math.min(hi, up[k] + h);
      const dn = u.slice(); dn[k] = Math.max(lo, dn[k] - h);
      grad[k] = (f(up) - f(dn)) / (up[k] - dn[k] || 2 * h);
    }
    const gnorm = Math.sqrt(grad.reduce((s, x) => s + x * x, 0)) || 1e-12;

    // backtracking Armijo sur la direction de descente projetée
    let accepted = false;
    let trialStep = step;
    for (let bt = 0; bt < 20; bt++) {
      const uTrial = projectBoxSum(
        u.map((x, k) => x - (trialStep / gnorm) * grad[k]), lo, hi, budget
      );
      const fTrial = f(uTrial);
      if (fTrial < fu - 1e-12) {
        u = uTrial; fu = fTrial; accepted = true;
        step = Math.min(trialStep * 1.3, opts.step0 || 0.5);
        break;
      }
      trialStep *= 0.5;
    }
    if (!accepted) {
      step *= 0.5;
      if (step < 1e-8) break;
    }
  }
  return { u, f: fu };
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = { optimiserAllocation, projectBoxSum };
}
