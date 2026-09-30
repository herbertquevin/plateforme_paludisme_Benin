// Port JS fidèle de scripts/{core_data,mobility,model}.py (mémoire "Modélisation
// de la propagation géographique du paludisme au Bénin"). Chaque fonction cite
// le(s) numéro(s) d'équation qu'elle implémente, comme dans le code Python source.
// Aucune formule numérique simplifiée, agrégée ou "équivalente" n'est utilisée :
// port direct, terme à terme, du modèle SEIRS-SEI multi-patchs à 12 départements.

(function (global) {
  "use strict";

  // -------------------------------------------------------------------
  // Données brutes par département (core_data.py)
  // -------------------------------------------------------------------
  const DEPARTEMENTS = [
    "Alibori", "Atacora", "Atlantique", "Borgou", "Collines", "Couffo",
    "Donga", "Littoral", "Mono", "Ouémé", "Plateau", "Zou",
  ];
  const N_PATCH = DEPARTEMENTS.length;

  const POPULATION = {
    Alibori: 867463, Atacora: 772262, Atlantique: 1398229, Borgou: 1214249,
    Collines: 717477, Couffo: 745328, Donga: 543130, Littoral: 679012,
    Mono: 497243, Ouémé: 1100404, Plateau: 622372, Zou: 851580,
  };

  const HBR = {
    Alibori: 11.08, Atacora: 15.09, Atlantique: 1.55, Borgou: 1.56,
    Collines: 19.41, Couffo: 2.78, Donga: 8.75, Littoral: 34.59,
    Mono: 6.91, Ouémé: 1.13, Plateau: 0.92, Zou: 9.22,
  };

  const PARTURITE = {
    Alibori: 81.0, Atacora: 85.0, Atlantique: 77.0, Borgou: 52.1,
    Collines: 77.8, Couffo: 90.0, Donga: 83.7, Littoral: 73.0,
    Mono: 73.4, Ouémé: 90.4, Plateau: 84.0, Zou: 81.0,
  };

  // Distances routières réelles entre chefs-lieux (km) -- matrice NOMINALE,
  // équation (2.7). Source : OpenStreetMap/OSRM (interrogé le 11 septembre
  // 2026), symétrisée par moyenne aller/retour. Issue de
  // sims/distances_routieres.json (ordre = DEPARTEMENTS).
  const _DISTANCES_ROUTIERES_BRUTES = [[0,252.0254,574.8756999999999,214.2287,426.04940000000005,558.3371,279.2201,629.658,590.192,574.7470999999999,502.6611,510.59009999999995],[251.9913,0,484.4043,215.289,335.57809999999995,467.8658,78.63260000000001,539.1866,499.7207,517.3249000000001,468.71590000000003,420.11879999999996],[574.844,484.3198,0,360.748,148.8956,93.0771,406.1088,55.1505,100.58810000000001,92.2239,169.351,82.3968],[214.4061,215.2547,360.7794,0,211.9531,344.2408,136.9295,415.56170000000003,376.0957,360.6508,288.5648,296.49379999999996],[426.01779999999997,335.49359999999996,148.9923,211.9218,0,132.4537,257.28270000000003,203.7746,164.30870000000002,181.9129,133.3039,84.7067],[560.0255999999999,469.5013,93.137,345.9296,134.0772,0,391.29040000000003,147.9193,36.856,184.9927,154.5325,46.3386],[279.2502,78.5988,406.3338,137.0043,257.5075,389.7952,0,461.11609999999996,421.65009999999995,439.2543,390.6454,342.0482],[630.0985,539.5741999999999,55.6049,416.0024,204.1501,148.3315,461.3633,0,108.1664,36.8339,97.7853,137.6513],[591.0663000000001,500.542,101.5398,376.9703,165.1179,37.9583,422.3311,108.5428,0,145.61620000000002,185.5733,77.3793],[575.0179,517.3311,91.4316,360.92190000000005,181.9069,184.15820000000002,439.1202,36.8337,143.9931,0,61.7764,173.478],[502.63190000000003,463.0564,163.856,288.5358,127.6322,148.97729999999999,384.8455,99.31230000000001,180.83229999999998,63.3664,0,104.95410000000001],[515.9811,425.4568,82.3328,301.88509999999997,90.03269999999999,46.296699999999994,347.2459,137.1151,78.1516,174.1885,110.488,0]]; // remplacé à l'inlining

  function _symmetrize(M) {
    const n = M.length;
    const out = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        out[i][j] = (M[i][j] + M[j][i]) / 2.0;
      }
    }
    return out;
  }

  const DISTANCES_KM = _symmetrize(_DISTANCES_ROUTIERES_BRUTES);

  // -------------------------------------------------------------------
  // Paramètres communs (params.py)
  // -------------------------------------------------------------------
  const MU_H = 4.29e-5;
  const SIGMA = 0.1;
  const NU = 0.0833;
  const A = 0.5;
  const BETA = 0.022;
  const B = 0.48;
  const ALPHA = 0.0043;
  const RHO = 0.0029;
  const OMEGA = 0.0005;
  const EPSILON = 0.2;
  const DELTA = 6.7e-6;
  const D_GONOTROPHIQUE = 2.5;
  const THETA_NOMINAL = 0.005;
  const KAPPA_NOM = 0.01479;

  // -------------------------------------------------------------------
  // Algèbre linéaire minimale (vecteurs = Array<number>, matrices = Array<Array<number>>)
  // -------------------------------------------------------------------
  const zeros = (n) => new Array(n).fill(0);
  const ones = (n) => new Array(n).fill(1);
  const add = (u, v) => u.map((x, i) => x + v[i]);
  const sub = (u, v) => u.map((x, i) => x - v[i]);
  const mul = (u, v) => u.map((x, i) => x * v[i]);
  const div = (u, v) => u.map((x, i) => x / v[i]);
  const scale = (u, s) => u.map((x) => x * s);
  const asArray = (x, n) => (Array.isArray(x) ? x : new Array(n).fill(x));

  // (P @ v)[k] = sum_j P[k][j] v[j]
  function matVec(P, v) {
    const n = P.length;
    const out = zeros(n);
    for (let k = 0; k < n; k++) {
      let s = 0;
      for (let j = 0; j < n; j++) s += P[k][j] * v[j];
      out[k] = s;
    }
    return out;
  }
  // (P^T @ v)[k] = sum_j P[j][k] v[j]
  function matVecT(P, v) {
    const n = P.length;
    const out = zeros(n);
    for (let k = 0; k < n; k++) {
      let s = 0;
      for (let j = 0; j < n; j++) s += P[j][k] * v[j];
      out[k] = s;
    }
    return out;
  }
  function matMul(A_, B_) {
    const n = A_.length, m = B_[0].length, p = B_.length;
    const out = Array.from({ length: n }, () => new Array(m).fill(0));
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < p; k++) {
        const aik = A_[i][k];
        if (aik === 0) continue;
        for (let j = 0; j < m; j++) out[i][j] += aik * B_[k][j];
      }
    }
    return out;
  }
  function matVecMul(M, v) {
    return matVec(M, v);
  }

  // -------------------------------------------------------------------
  // Matrice de mobilité (mobility.py) -- (2.7)-(2.8)
  // -------------------------------------------------------------------
  function matriceResidence(theta = THETA_NOMINAL, expDist = 2.0) {
    const Nh = DEPARTEMENTS.map((d) => POPULATION[d]);
    const n = N_PATCH;
    const d = DISTANCES_KM.map((row) => row.slice());
    for (let i = 0; i < n; i++) d[i][i] = 1.0; // évite division par zéro (annulé juste après)

    const W = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) {
        W[k][j] = (Nh[k] * Nh[j]) / Math.pow(d[k][j], expDist); // (2.7)
      }
    }
    for (let i = 0; i < n; i++) W[i][i] = 0.0;

    const colSums = new Array(n).fill(0);
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += W[k][j];
      colSums[j] = s;
    }
    const M = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) M[k][j] = W[k][j] / colSums[j]; // (2.8)
    }

    const P = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) {
        P[k][j] = (1.0 - theta) * (k === j ? 1 : 0) + theta * M[k][j]; // (2.8)
      }
    }
    return { M, P };
  }

  // -------------------------------------------------------------------
  // Paramétrisation vectorielle départementale (model.py)
  // -------------------------------------------------------------------
  function muVVector(D = D_GONOTROPHIQUE) {
    // mu_{v,k} = -ln(parturite_k) / D (Detinova, 1962)
    return DEPARTEMENTS.map((d_) => -Math.log(PARTURITE[d_] / 100.0) / D);
  }

  function NVVector(kappa, muV, a = A) {
    const Nh = DEPARTEMENTS.map((d_) => POPULATION[d_]);
    const hbr = DEPARTEMENTS.map((d_) => HBR[d_]);
    return Nh.map((nh, k) => kappa * (hbr[k] / a) * nh);
  }

  function buildParams(kappa = KAPPA_NOM, overrides = {}) {
    const a = overrides.a !== undefined ? overrides.a : A;
    const beta = overrides.beta !== undefined ? overrides.beta : BETA;
    const b = overrides.b !== undefined ? overrides.b : B;
    const sigma_ = overrides.sigma !== undefined ? overrides.sigma : SIGMA;
    const nu_ = overrides.nu !== undefined ? overrides.nu : NU;
    const alpha_ = overrides.alpha !== undefined ? overrides.alpha : ALPHA;
    const rho_ = overrides.rho !== undefined ? overrides.rho : RHO;
    const omega_ = overrides.omega !== undefined ? overrides.omega : OMEGA;
    const epsilon_ = overrides.epsilon !== undefined ? overrides.epsilon : EPSILON;
    const delta_ = overrides.delta !== undefined ? overrides.delta : DELTA;
    const mu_ = overrides.mu !== undefined ? overrides.mu : MU_H;

    const N_h = DEPARTEMENTS.map((d_) => POPULATION[d_]);
    const mu_v = muVVector();
    const N_v = NVVector(kappa, mu_v, a);
    const Lambda_v = mul(mu_v, N_v); // (2.18)

    const n = N_PATCH;
    return {
      N_h, N_v, mu_v, Lambda_v,
      a, beta, b,
      sigma: ones(n).map(() => sigma_),
      nu: ones(n).map(() => nu_),
      alpha: ones(n).map(() => alpha_),
      rho: ones(n).map(() => rho_),
      omega: ones(n).map(() => omega_),
      epsilon: epsilon_,
      delta: ones(n).map(() => delta_),
      mu: ones(n).map(() => mu_),
      Lambda_h: N_h.map((nh) => mu_ * nh), // (2.18')
    };
  }

  // -------------------------------------------------------------------
  // Système différentiel complet (2.9), forces d'infection (2.4)-(2.5),
  // présence humaine effective (2.2)
  // -------------------------------------------------------------------
  function rhsMultipatch(y, p, P) {
    const n = N_PATCH;
    const Sh = y.slice(0, n), Eh = y.slice(n, 2 * n), Ih = y.slice(2 * n, 3 * n), Rh = y.slice(3 * n, 4 * n);
    const Sv = y.slice(4 * n, 5 * n), Ev = y.slice(5 * n, 6 * n), Iv = y.slice(6 * n, 7 * n);

    const N_h = add(add(Sh, Eh), add(Ih, Rh)); // (2.1)
    const N_h_tilde = matVec(P, N_h); // (2.2)
    const N_h_tilde_safe = N_h_tilde.map((x) => Math.max(x, 1e-12));

    const aArr = asArray(p.a, n);
    // (2.4) lambda_h^k = (P^T @ (a*beta*Iv/Ñh))_k
    const vTerm = Iv.map((iv, j) => (aArr[j] * p.beta * iv) / N_h_tilde_safe[j]);
    const lambda_h = matVecT(P, vTerm);

    // (2.3) I_tilde_h^k = sum_j p_kj (Ih^j + eps Rh^j)
    const IepsR = Ih.map((ih, j) => ih + p.epsilon * Rh[j]);
    const I_tilde_h = matVec(P, IepsR);
    // (2.5) lambda_v^k = a*b*I_tilde_h^k / Ñh^k
    const lambda_v = I_tilde_h.map((it, k) => (aArr[k] * p.b * it) / N_h_tilde_safe[k]);

    const dSh = zeros(n), dEh = zeros(n), dIh = zeros(n), dRh = zeros(n);
    const dSv = zeros(n), dEv = zeros(n), dIv = zeros(n);
    for (let k = 0; k < n; k++) {
      dSh[k] = p.Lambda_h[k] + p.rho[k] * Ih[k] + p.omega[k] * Rh[k] - (lambda_h[k] + p.mu[k]) * Sh[k];
      dEh[k] = lambda_h[k] * Sh[k] - (p.sigma[k] + p.mu[k]) * Eh[k];
      dIh[k] = p.sigma[k] * Eh[k] - (p.alpha[k] + p.rho[k] + p.delta[k] + p.mu[k]) * Ih[k];
      dRh[k] = p.alpha[k] * Ih[k] - (p.omega[k] + p.mu[k]) * Rh[k];
      dSv[k] = p.Lambda_v[k] - (lambda_v[k] + p.mu_v[k]) * Sv[k];
      dEv[k] = lambda_v[k] * Sv[k] - (p.nu[k] + p.mu_v[k]) * Ev[k];
      dIv[k] = p.nu[k] * Ev[k] - p.mu_v[k] * Iv[k];
    }
    return [...dSh, ...dEh, ...dIh, ...dRh, ...dSv, ...dEv, ...dIv];
  }

  // -------------------------------------------------------------------
  // Matrice de prochaine génération spatiale (2.22)-(2.23)
  // -------------------------------------------------------------------
  function KhvMatrix(p, P) {
    const n = N_PATCH;
    const aArr = asArray(p.a, n);
    const N_h0 = p.N_h, N_v0 = p.N_v, mu_v = p.mu_v, nu = p.nu;
    const sigma = p.sigma, mu = p.mu, alpha = p.alpha, rho_ = p.rho, delta = p.delta, omega = p.omega;
    const epsilon = p.epsilon, beta = p.beta, b = p.b;

    const N_h_tilde0 = matVec(P, N_h0); // (2.2) à l'ESM

    const A_ell = new Array(n);
    for (let ell = 0; ell < n; ell++) {
      A_ell[ell] = (aArr[ell] * aArr[ell] * b * beta * N_v0[ell] * nu[ell]) /
        (mu_v[ell] * (nu[ell] + mu_v[ell]) * N_h_tilde0[ell] * N_h_tilde0[ell]);
    }
    const B_j = new Array(n);
    for (let j = 0; j < n; j++) {
      B_j[j] = (sigma[j] / ((sigma[j] + mu[j]) * (alpha[j] + rho_[j] + delta[j] + mu[j]))) *
        (1.0 + (epsilon * alpha[j]) / (omega[j] + mu[j]));
    }

    // Q[k][j] = sum_ell P[ell][k] * A_ell[ell] * P[ell][j]  (P^T diag(A) P)
    const Q = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let ell = 0; ell < n; ell++) {
      const a_ell = A_ell[ell];
      if (a_ell === 0) continue;
      for (let k = 0; k < n; k++) {
        const pek = P[ell][k];
        if (pek === 0) continue;
        const coeff = a_ell * pek;
        for (let j = 0; j < n; j++) Q[k][j] += coeff * P[ell][j];
      }
    }

    const K_hv = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let k = 0; k < n; k++) {
      for (let j = 0; j < n; j++) K_hv[k][j] = N_h0[k] * Q[k][j] * B_j[j];
    }
    return K_hv;
  }

  // Rayon spectral d'une matrice non négative par itération de la puissance
  // (Perron-Frobenius : valeur propre dominante réelle et positive).
  function spectralRadiusNonNeg(M, iters = 400, tol = 1e-10) {
    const n = M.length;
    let v = ones(n).map((x) => x / n);
    let lambda = 0;
    for (let it = 0; it < iters; it++) {
      const w = matVec(M, v);
      const norm = Math.sqrt(w.reduce((s, x) => s + x * x, 0));
      if (norm < 1e-300) return 0;
      const vNext = w.map((x) => x / norm);
      const lambdaNext = norm; // ~ ||M v|| quand v normalisé (converge vers rho)
      if (Math.abs(lambdaNext - lambda) < tol * Math.max(1, Math.abs(lambdaNext))) {
        v = vNext;
        lambda = lambdaNext;
        break;
      }
      v = vNext;
      lambda = lambdaNext;
    }
    // Rayleigh quotient final pour affiner : lambda = (v^T M v) / (v^T v)
    const Mv = matVec(M, v);
    const num = v.reduce((s, x, i) => s + x * Mv[i], 0);
    const den = v.reduce((s, x) => s + x * x, 0);
    return den > 0 ? num / den : lambda;
  }

  function R0Spatial(p, P) {
    const K_hv = KhvMatrix(p, P);
    const rho = spectralRadiusNonNeg(K_hv);
    return Math.sqrt(Math.max(rho, 0.0));
  }

  function identity(n) {
    return Array.from({ length: n }, (_, i) => {
      const row = new Array(n).fill(0);
      row[i] = 1;
      return row;
    });
  }

  function R0LocalNoMobility(p) {
    const n = N_PATCH;
    const P_id = identity(n);
    const K_hv = KhvMatrix(p, P_id);
    return K_hv.map((row, k) => Math.sqrt(Math.max(row[k], 0.0)));
  }

  function R0DiagAvecMobilite(p, P) {
    const K_hv = KhvMatrix(p, P);
    return K_hv.map((row, k) => Math.sqrt(Math.max(row[k], 0.0)));
  }

  function sourcePuitsIndex(p, P) {
    const K_hv = KhvMatrix(p, P);
    const n = K_hv.length;
    const I_k = K_hv.map((row) => row.reduce((s, x) => s + x, 0)); // sum axis=1 (importation)
    const E_k = new Array(n).fill(0);
    for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) E_k[j] += K_hv[k][j]; // sum axis=0 (exportation)
    return E_k.map((e, k) => e - I_k[k]);
  }

  // -------------------------------------------------------------------
  // Scénarios de contrôle (S11-S15)
  // -------------------------------------------------------------------
  function controleA(p, u) {
    const n = N_PATCH;
    const uArr = asArray(u, n);
    const aArr = asArray(p.a, n);
    return Object.assign({}, p, { a: aArr.map((av, k) => (1.0 - uArr[k]) * av) });
  }

  function controleLambdaV(p, u) {
    const n = N_PATCH;
    const uArr = asArray(u, n);
    return Object.assign({}, p, {
      Lambda_v: p.Lambda_v.map((lv, k) => (1.0 - uArr[k]) * lv),
      N_v: p.N_v.map((nv, k) => (1.0 - uArr[k]) * nv),
    });
  }

  // -------------------------------------------------------------------
  // Conditions initiales et intégration (RK4 à pas fixe jusqu'à l'équilibre)
  // -------------------------------------------------------------------
  function etatInitialPrevalence(p, prevH, fracE = 0.10, prevV = 0.002, fracEv = 0.001) {
    const n = N_PATCH;
    const N_h = p.N_h, N_v = p.N_v;
    const Ih0 = N_h.map((nh) => prevH * nh);
    const Eh0 = Ih0.map((ih) => fracE * ih);
    const Rh0 = zeros(n);
    const Sh0 = N_h.map((nh, k) => nh - Eh0[k] - Ih0[k] - Rh0[k]);
    const Iv0 = N_v.map((nv) => prevV * nv);
    const Ev0 = N_v.map((nv) => fracEv * nv);
    const Sv0 = N_v.map((nv, k) => nv - Ev0[k] - Iv0[k]);
    return [...Sh0, ...Eh0, ...Ih0, ...Rh0, ...Sv0, ...Ev0, ...Iv0];
  }

  function rk4Step(y, p, P, dt) {
    const addScaled = (u, v, s) => u.map((x, i) => x + s * v[i]);
    const k1 = rhsMultipatch(y, p, P);
    const k2 = rhsMultipatch(addScaled(y, k1, dt / 2), p, P);
    const k3 = rhsMultipatch(addScaled(y, k2, dt / 2), p, P);
    const k4 = rhsMultipatch(addScaled(y, k3, dt), p, P);
    const out = new Array(y.length);
    for (let i = 0; i < y.length; i++) {
      out[i] = y[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    return out;
  }

  // Intègre (2.9) jusqu'à t_end_years par RK4 à pas fixe (dt en jours).
  function simuler(y0, p, P, tEndYears = 150.0, dt = 1.0) {
    const tEnd = 365.25 * tEndYears;
    let y = y0.slice();
    let t = 0;
    const nSteps = Math.round(tEnd / dt);
    for (let i = 0; i < nSteps; i++) {
      y = rk4Step(y, p, P, dt);
      t += dt;
    }
    return { y, t };
  }

  function extraire(y) {
    const n = N_PATCH;
    const Sh = y.slice(0, n), Eh = y.slice(n, 2 * n), Ih = y.slice(2 * n, 3 * n), Rh = y.slice(3 * n, 4 * n);
    const Sv = y.slice(4 * n, 5 * n), Ev = y.slice(5 * n, 6 * n), Iv = y.slice(6 * n, 7 * n);
    const N_h = add(add(Sh, Eh), add(Ih, Rh));
    const N_v = add(add(Sv, Ev), Iv);
    return { Sh, Eh, Ih, Rh, Sv, Ev, Iv, N_h, N_v };
  }

  // Résultats d'une analyse Monte-Carlo (2000 tirages) sur les incertitudes des
  // paramètres entomologiques/immunitaires — voir onglet Robustesse et carte
  // "Fiabilité du classement" de la page Résultats clés.
  const ROBUSTESSE_MC = {
    Atacora:{pr_source:100.0,pr_top3:96.5,r0_local:1.9046}, Littoral:{pr_source:100.0,pr_top3:92.5,r0_local:1.7448},
    Collines:{pr_source:99.9,pr_top3:64.5,r0_local:1.5616}, Donga:{pr_source:94.8,pr_top3:19.6,r0_local:1.3582},
    Alibori:{pr_source:94.8,pr_top3:16.0,r0_local:1.3483}, Zou:{pr_source:87.3,pr_top3:3.8,r0_local:1.2299},
    Couffo:{pr_source:66.8,pr_top3:7.2,r0_local:1.1040}, Ouémé:{pr_source:10.4,pr_top3:0.0,r0_local:0.7243},
    Mono:{pr_source:5.3,pr_top3:0.0,r0_local:0.7909}, Atlantique:{pr_source:0.0,pr_top3:0.0,r0_local:0.4277},
    Borgou:{pr_source:0.0,pr_top3:0.0,r0_local:0.2007}, Plateau:{pr_source:0.0,pr_top3:0.0,r0_local:0.4470},
  };
  const N_TIRAGES_MC = 2000;

  const MalariaModel = {
    DEPARTEMENTS, N_PATCH, POPULATION, HBR, PARTURITE, DISTANCES_KM,
    MU_H, SIGMA, NU, A, BETA, B, ALPHA, RHO, OMEGA, EPSILON, DELTA,
    D_GONOTROPHIQUE, THETA_NOMINAL, KAPPA_NOM,
    ROBUSTESSE_MC, N_TIRAGES_MC,
    matriceResidence, buildParams, rhsMultipatch,
    KhvMatrix, R0Spatial, R0LocalNoMobility, R0DiagAvecMobilite, sourcePuitsIndex,
    controleA, controleLambdaV,
    etatInitialPrevalence, simuler, extraire,
    _linalg: { matVec, matVecT, matMul, add, sub, mul, div, scale, asArray, spectralRadiusNonNeg },
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = MalariaModel;
  } else {
    global.MalariaModel = MalariaModel;
  }
})(typeof window !== "undefined" ? window : globalThis);
