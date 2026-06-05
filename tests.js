// tests.js — pure Node tests for engine.js
// Every expected value is hand-derived independently (NOT by running the engine).
// Run: node tests.js   ->  exits 0 if all pass.

const { Engine, Metric } = require('./engine.js');

let N = 0, M = 0;
function approx(label, got, exp, tol) {
    tol = (tol === undefined) ? 1e-6 : tol;
    if (typeof got === 'number' && Math.abs(got - exp) <= tol) { N++; console.log('PASS ' + label + '  (' + got + ')'); }
    else { M++; console.log('FAIL ' + label + '  got=' + got + ' expected=' + exp + ' tol=' + tol); }
}
function ok(label, cond, info) {
    if (cond) { N++; console.log('PASS ' + label); }
    else { M++; console.log('FAIL ' + label + (info ? '  ' + info : '')); }
}

// ---------------------------------------------------------------------------
// 1. Normal-tail p (two-sided). getP is a p-value, NOT a CDF.
//    At z=0 the two-sided tail probability is 1.0 (NOT 0.5, NOT 0).
// ---------------------------------------------------------------------------
approx('getP(0) two-sided p = 1.0', Engine.getP(0), 1.0, 1e-3);
//    At z=1.96 the two-sided p is ~0.05.
approx('getP(1.96) ~= 0.05', Engine.getP(1.96), 0.05, 5e-4);
//    Symmetry: getP uses |z|.
approx('getP(-1.96) == getP(1.96)', Engine.getP(-1.96), Engine.getP(1.96), 1e-12);

// ---------------------------------------------------------------------------
// 2. t-approximation getT(df) = 1.96 + 2.37/df  (df>0); 1.96 for df<=0.
//    Hand: getT(1)=1.96+2.37=4.33 ; getT(10)=1.96+0.237=2.197 ; getT(0)=1.96.
// ---------------------------------------------------------------------------
approx('getT(1) = 4.33', Engine.getT(1), 4.33, 1e-9);
approx('getT(10) = 2.197', Engine.getT(10), 2.197, 1e-9);
approx('getT(0) = 1.96 (floor)', Engine.getT(0), 1.96, 1e-12);

// ---------------------------------------------------------------------------
// 3. FULLY HAND-WORKED 2-STUDY POOLING (Continuous, MD, fixed, IV).
//    Study A: m1=10,s1=2,n1=10 ; m2=5,s2=2,n2=10
//        MD_A = 10-5 = 5 ;  vi_A = 2^2/10 + 2^2/10 = 0.4+0.4 = 0.8 ;  w_A = 1/0.8 = 1.25
//    Study B: m1=12,s1=3,n1=20 ; m2=6,s2=3,n2=20
//        MD_B = 12-6 = 6 ;  vi_B = 9/20 + 9/20 = 0.45+0.45 = 0.9 ;  w_B = 1/0.9 = 1.111111...
//    Pooled (fixed): mu = (1.25*5 + 1.111111*6)/(1.25+1.111111)
//                       = (6.25 + 6.666667)/2.361111 = 12.916667/2.361111 = 5.4705882
//    pse = sqrt(1/2.361111) = sqrt(0.4235294) = 0.6507914
//    z = 1.96 (HKSJ off) -> lo = 5.4705882 - 1.96*0.6507914 = 4.1950371
//                            hi = 5.4705882 + 1.96*0.6507914 = 6.7461393
//    Q = 1.25*(5-5.4705882)^2 + 1.111111*(6-5.4705882)^2
//      = 1.25*0.2214533 + 1.111111*0.2802768 = 0.2768166 + 0.3114187 = 0.5882353
//    df = 1 ; tau^2(DL,fixed) = 0 ; I^2 = max(0,(0.5882353-1)/0.5882353)*100 = 0
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10, s1: 2, n1: 10, m2: 5, s2: 2, n2: 10 },
         { id: 'B', m1: 12, s1: 3, n1: 20, m2: 6, s2: 3, n2: 20 }], 'cont', 'MD');
    approx('MD A es', eff[0].es, 5, 1e-12);
    approx('MD A vi', eff[0].vi, 0.8, 1e-12);
    approx('MD B es', eff[1].es, 6, 1e-12);
    approx('MD B vi', eff[1].vi, 0.9, 1e-12);
    const r = Engine.pool(eff, 'fixed', 'iv', 'dl', false);
    approx('MD pooled es', r.es, 5.4705882353, 1e-9);
    approx('MD pooled lo', r.lo, 4.1950371433, 1e-7);
    approx('MD pooled hi', r.hi, 6.7461393273, 1e-7);
    approx('MD pooled displayVal == es (MD is linear)', r.displayVal, r.es, 1e-12);
    approx('MD Q', r.Q, 0.5882353, 1e-6);
    approx('MD tau2 (fixed=0)', r.tau2, 0, 1e-12);
    approx('MD I2 = 0', r.I2, 0, 1e-9);
}

// ---------------------------------------------------------------------------
// 4. EDGE: two identical studies => tau^2 = 0, I^2 = 0, Q = 0.
//    OR, e1=10,n1=100,e2=20,n2=100 (no zero cells => cc=0).
//    OR = (10*80)/(20*90) = 800/1800 = 0.444444 ; logOR = ln(0.444444) = -0.8109302
//    vi = 1/10 + 1/90 + 1/20 + 1/80 = 0.1+0.0111111+0.05+0.0125 = 0.1736111
//    Two identical => Q=0, tau^2=0, I^2=0 ; pooled logOR = -0.8109302 ; disp=exp(.)=0.444444
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: '1', e1: 10, n1: 100, e2: 20, n2: 100 },
         { id: '2', e1: 10, n1: 100, e2: 20, n2: 100 }], 'binary', 'OR');
    approx('OR each es (logOR)', eff[0].es, -0.8109302162, 1e-9);
    approx('OR each vi', eff[0].vi, 0.1736111111, 1e-9);
    const r = Engine.pool(eff, 'random', 'iv', 'dl', false);
    approx('OR identical Q = 0', r.Q, 0, 1e-12);
    approx('OR identical tau2 = 0', r.tau2, 0, 1e-12);
    approx('OR identical I2 = 0', r.I2, 0, 1e-12);
    approx('OR pooled logOR', r.es, -0.8109302162, 1e-9);
    approx('OR pooled displayVal (back-transform)', r.displayVal, 0.4444444444, 1e-9);
    ok('OR pooled on LOG scale then back-transformed', Math.abs(Math.exp(r.es) - r.displayVal) < 1e-12);
}

// ---------------------------------------------------------------------------
// 5. EDGE: k=1 passthrough.
//    Single MD study: m1=10,s1=2,n1=10,m2=5,s2=2,n2=10 => es=5.
//    Q = 0, df = 0, I^2 = (0-0)/0 = NaN (0/0). Pooled es passes through = 5.
//    No prediction interval (k<=2) => pi_lo null.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'X', m1: 10, s1: 2, n1: 10, m2: 5, s2: 2, n2: 10 }], 'cont', 'MD');
    ok('k=1 produces one effect', eff.length === 1);
    const r = Engine.pool(eff, 'random', 'iv', 'dl', true);
    approx('k=1 pooled es passthrough', r.es, 5, 1e-12);
    approx('k=1 displayVal passthrough', r.displayVal, 5, 1e-12);
    ok('k=1 I2 is NaN (0/0)', Number.isNaN(r.I2));
    ok('k=1 no prediction interval', r.pi_lo === null && r.pi_hi === null);
}

// ---------------------------------------------------------------------------
// 6. EDGE: empty guards.
//    calcEffects([]) -> [] ; pool([]) -> null.
// ---------------------------------------------------------------------------
ok('calcEffects([]) returns empty array', Array.isArray(Engine.calcEffects([], 'binary', 'OR')) && Engine.calcEffects([], 'binary', 'OR').length === 0);
ok('pool([]) returns null', Engine.pool([], 'fixed', 'iv', 'dl', false) === null);
ok('pool(all-excluded) returns null', Engine.pool([{ es: 1, vi: 1, excluded: true, metric: 'MD' }], 'fixed', 'iv', 'dl', false) === null);

// ---------------------------------------------------------------------------
// 7. HKSJ floor + t-critical (HAND-WORKED, 2 studies, Continuous MD, random).
//    A: MD=5, vi=4/30+4/30 = 0.2666667 ; B: MD=1, vi=0.2666667.
//    Fixed weights w=1/0.2666667=3.75 ; sw=7.5 ; mu=(3.75*5+3.75*1)/7.5 = 22.5/7.5 = 3.
//    Q = 3.75*(5-3)^2 + 3.75*(1-3)^2 = 15+15 = 30 ; df=1.
//    C = sw - sum(w^2)/sw = 7.5 - (3.75^2*2)/7.5 = 7.5 - 28.125/7.5 = 7.5-3.75 = 3.75.
//    tau^2(DL) = (30-1)/3.75 = 29/3.75 = 7.7333333.
//    RE weights w* = 1/(0.2666667+7.7333333) = 1/8 = 0.125 ; s = 0.25.
//    pes = (0.125*5 + 0.125*1)/0.25 = 0.75/0.25 = 3.
//    pse_RE = sqrt(1/0.25) = 2.
//    vhk = [0.125*(5-3)^2 + 0.125*(1-3)^2] / ((k-1)*s) = [0.5+0.5]/(1*0.25) = 4.
//    floor = max(pse_RE^2=4, vhk=4) = 4 (q=Q/df=1 here, floor binds at equality) -> pse=2.
//    z = getT(df=1) = 1.96 + 2.37 = 4.33  (t_{k-1} critical via approximation, NOT 1.96).
//    lo = 3 - 4.33*2 = -5.66 ; hi = 3 + 4.33*2 = 11.66.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10, s1: 2, n1: 30, m2: 5, s2: 2, n2: 30 },
         { id: 'B', m1: 8, s1: 2, n1: 30, m2: 7, s2: 2, n2: 30 }], 'cont', 'MD');
    approx('HKSJ A vi', eff[0].vi, 0.2666666667, 1e-9);
    const r = Engine.pool(eff, 'random', 'iv', 'dl', true);
    approx('HKSJ tau2', r.tau2, 7.7333333333, 1e-7);
    approx('HKSJ pooled es', r.es, 3, 1e-9);
    approx('HKSJ se (with floor)', r.se, 2, 1e-9);
    approx('HKSJ lo (t-critical 4.33)', r.lo, -5.66, 1e-7);
    approx('HKSJ hi (t-critical 4.33)', r.hi, 11.66, 1e-7);
    // Confirm HKSJ used a t-critical, NOT z=1.96: a z-based CI would be 3 +/- 1.96*2 = [-0.92, 6.92].
    ok('HKSJ widened beyond z (not z=1.96 in HKSJ)', r.lo < -0.92 - 1e-6);
}

// ---------------------------------------------------------------------------
// 8. HKSJ variance FLOOR: when q = Q/df < 1, the floor must keep RE variance.
//    Two nearly-identical studies => Q < df => HKSJ se must NOT shrink below RE se.
//    A,B almost identical (MD 5.0 vs 5.0) but tiny diff -> q<1 -> floor binds at RE se.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10.00, s1: 2, n1: 30, m2: 5, s2: 2, n2: 30 },
         { id: 'B', m1: 10.01, s1: 2, n1: 30, m2: 5, s2: 2, n2: 30 }], 'cont', 'MD');
    const rRE = Engine.pool(eff, 'random', 'iv', 'dl', false);
    const rHK = Engine.pool(eff, 'random', 'iv', 'dl', true);
    // RE se with z=1.96; HKSJ se must be >= RE se (floor), never below.
    ok('HKSJ se floored at >= RE se (no shrink below 1)', rHK.se >= rRE.se - 1e-12,
        'rHK.se=' + rHK.se + ' rRE.se=' + rRE.se);
}

// ---------------------------------------------------------------------------
// 9. BCG vaccine benchmark (RR, random, DL). Classic dataset.
//    Engine uses normal-approx log-RR; tolerance 0.02 vs metafor reference.
//    Reference (metafor): pooled logRR ~= -0.7141, tau^2 ~= 0.3088.
// ---------------------------------------------------------------------------
{
    const bcg = [{ id: '1', e1: 4, n1: 123, e2: 11, n2: 139 }, { id: '2', e1: 6, n1: 306, e2: 29, n2: 303 }, { id: '3', e1: 3, n1: 231, e2: 11, n2: 220 }, { id: '4', e1: 62, n1: 13598, e2: 248, n2: 12867 }, { id: '5', e1: 33, n1: 5069, e2: 47, n2: 5808 }, { id: '6', e1: 180, n1: 1541, e2: 372, n2: 1451 }, { id: '7', e1: 8, n1: 2545, e2: 10, n2: 629 }, { id: '8', e1: 505, n1: 88391, e2: 499, n2: 88391 }, { id: '9', e1: 29, n1: 7499, e2: 45, n2: 7277 }, { id: '10', e1: 17, n1: 1716, e2: 65, n2: 1617 }, { id: '11', e1: 186, n1: 50634, e2: 141, n2: 27338 }, { id: '12', e1: 5, n1: 2498, e2: 3, n2: 2341 }, { id: '13', e1: 27, n1: 16913, e2: 29, n2: 17854 }];
    const eff = Engine.calcEffects(bcg, 'binary', 'RR');
    const r = Engine.pool(eff, 'random', 'iv', 'dl', false);
    approx('BCG pooled logRR ~ -0.714', r.es, -0.7141, 0.02);
    approx('BCG tau2 ~ 0.309', r.tau2, 0.3088, 0.02);
    ok('BCG pooled RR (back-transformed) < 1', r.displayVal < 1 && Math.abs(Math.exp(r.es) - r.displayVal) < 1e-9);
}

// ---------------------------------------------------------------------------
// 10. SMD (continuous) uses Hedges' g small-sample correction.
//     A: m1=10,s1=2,n1=20 ; m2=8,s2=2,n2=20. df=38.
//     sp = sqrt(((19*4)+(19*4))/38) = sqrt(152/38) = sqrt(4) = 2.
//     d  = (10-8)/2 = 1.0 ; J = 1 - 3/(4*38-1) = 1 - 3/151 = 0.9801325.
//     g  = 1.0 * 0.9801325 = 0.9801325.
//     vi = (n1+n2)/(n1*n2) + g^2/(2*(n1+n2)) = 40/400 + 0.9606597/80 = 0.1 + 0.0120082 = 0.1120082.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects(
        [{ id: 'A', m1: 10, s1: 2, n1: 20, m2: 8, s2: 2, n2: 20 }], 'cont', 'SMD');
    approx('SMD Hedges g', eff[0].es, 0.9801324503, 1e-9);
    approx('SMD vi', eff[0].vi, 0.1120082, 1e-6);
}

// ---------------------------------------------------------------------------
// 11. Survival HR: log scale, se from CI width / 3.92.
//     hr=0.8, lci=0.6, uci=0.95 => es=ln(0.8)=-0.2231436.
//     se = (ln(0.95)-ln(0.6))/3.92 = (-0.0512933 - (-0.5108256))/3.92 = 0.4595324/3.92 = 0.1172276.
//     vi = se^2 = 0.0137423.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects([{ id: '1', hr: 0.8, lci: 0.6, uci: 0.95, n1: 100 }], 'survival', 'HR');
    approx('HR es = ln(0.8)', eff[0].es, -0.2231435513, 1e-9);
    approx('HR vi', eff[0].vi, 0.0137423, 1e-6);
}

// ---------------------------------------------------------------------------
// 12. Proportion (logit). e=30,n=100 (no CC). p=0.3.
//     es = ln(0.3/0.7) = ln(0.4285714) = -0.8472979.
//     vi = 1/(n*p*(1-p)) = 1/(100*0.3*0.7) = 1/21 = 0.0476190.
//     displayVal = e/n = 0.3.
// ---------------------------------------------------------------------------
{
    const eff = Engine.calcEffects([{ id: '1', e: 30, n: 100 }], 'prop', 'Prop');
    approx('Prop logit es', eff[0].es, -0.8472978604, 1e-9);
    approx('Prop vi', eff[0].vi, 0.0476190476, 1e-9);
    approx('Prop displayVal', eff[0].displayVal, 0.3, 1e-12);
}

console.log('\n' + N + ' passed, ' + M + ' failed');
process.exit(M === 0 ? 0 : 1);
