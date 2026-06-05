// engine.js — pure statistical core for 786-MIII pairwise meta-analysis
// Extracted VERBATIM from the inline app script (single source of truth).
// Pure functions only: effect-size conversion, log-scale pooling (DL/REML/PM/SJ),
// tau^2 / I^2 / Q, CI, prediction interval, normal-tail p, t-approximation, Egger's test.
// No DOM / Plotly dependencies.

        const Metric={isLogRatio:m=>['RR','OR','HR'].includes(m)};
        const Engine={
            toNum:v=>{const x=parseFloat(v);return Number.isFinite(x)?x:0;},
            validNum:v=>Number.isFinite(v)&&!isNaN(v),
            getT:df=>(df<=0)?1.96:1.96+2.37/df,
            getP:z=>{z=Math.abs(z);const t=1/(1+0.2316419*z);return 2*(0.3989423*Math.exp(-0.5*z*z)*((((1.330274*t-1.821256)*t+1.781478)*t-0.356564)*t+0.319382)*t);},
            solveTau2:(eff,type,mod)=>{
                if(mod==='fixed')return 0;
                let sumW=0,sumWY=0; eff.forEach(e=>{const w=1/e.vi;sumW+=w;sumWY+=w*e.es;});
                const mu=sumWY/sumW, Q=eff.reduce((a,b)=>a+(1/b.vi)*(b.es-mu)**2,0), df=Math.max(0,eff.length-1);
                if(type==='dl'){ const swsq=eff.reduce((a,b)=>a+(1/b.vi)**2,0); return Math.max(0,(Q-df)/(sumW-swsq/sumW||1)); }
                if(type==='sj') return Math.max(0,eff.reduce((a,b)=>a+(b.es-mu)**2,0)/(df||1));
                let t2=0.01; for(let i=0;i<50;i++){
                    let w=eff.map(e=>1/(e.vi+t2)), sw=w.reduce((a,b)=>a+b,0), sw2=w.reduce((a,b)=>a+b*b,0);
                    let m=w.reduce((a,b,idx)=>a+b*eff[idx].es,0)/sw;
                    let step=0;
                    if(type==='reml'){ let qr=w.reduce((a,b,idx)=>a+b*(eff[idx].es-m)**2,0); step=(qr-(eff.length-sw2/sw))/w.reduce((a,b)=>a+b*b,0); }
                    else { let est=w.reduce((a,b,idx)=>a+b*(eff[idx].es-m)**2,0)-df; step=est/w.reduce((a,b,idx)=>a+b*b*(eff[idx].es-m)**2,0); }
                    if(Math.abs(step)>1)step=Math.sign(step); t2=Math.max(0,t2+step); if(Math.abs(step)<1e-6)break;
                }
                return t2;
            },
            calcEffects:(rows,type,met)=>{
                return rows.map((r,i)=>{
                    let es=0,vi=0,dv=0,lo=0,hi=0;
                    const n1=Engine.toNum(r.n1), n2=Engine.toNum(r.n2);
                    if(type==='binary'){
                        const e1=Engine.toNum(r.e1), e2=Engine.toNum(r.e2), cc=(e1===0||e2===0||e1===n1||e2===n2)?.5:0;
                        if(!n1||!n2)return null;
                        if(met==='RR'){ es=Math.log(((e1+cc)/(n1+cc))/((e2+cc)/(n2+cc))); vi=1/(e1+cc)-1/(n1+cc)+1/(e2+cc)-1/(n2+cc); dv=Math.exp(es); }
                        else if(met==='OR'){ es=Math.log(((e1+cc)*(n2-e2+cc))/((e2+cc)*(n1-e1+cc))); vi=1/(e1+cc)+1/(n1-e1+cc)+1/(e2+cc)+1/(n2-e2+cc); dv=Math.exp(es); }
                        else { const p1=(e1+cc)/(n1+2*cc), p2=(e2+cc)/(n2+2*cc); es=p1-p2; vi=(p1*(1-p1))/n1+(p2*(1-p2))/n2; dv=es; }
                    } else if(type==='cont'){
                        const m1=Engine.toNum(r.m1), s1=Engine.toNum(r.s1), m2=Engine.toNum(r.m2), s2=Engine.toNum(r.s2);
                        if(!n1||!n2)return null;
                        if(met==='MD'){ es=m1-m2; vi=(s1*s1)/n1+(s2*s2)/n2; dv=es; }
                        else { const df=n1+n2-2, sp=Math.sqrt(((n1-1)*s1*s1+(n2-1)*s2*s2)/df); es=(m1-m2)/sp*(1-(3/(4*df-1))); vi=(n1+n2)/(n1*n2)+(es*es)/(2*(n1+n2)); dv=es; }
                    } else if(type==='prop'){
                        const n=Engine.toNum(r.n), e=Engine.toNum(r.e); if(!n)return null;
                        const useCC = (e===0 || e===n);
                        const p_raw = e/n;
                        if(useCC){ const p_cc = (e+0.5)/(n+1); es = Math.log(p_cc/(1-p_cc)); vi = 1/((e+0.5)*(1-(e+0.5)/(n+1))); }
                        else { es = Math.log(p_raw/(1-p_raw)); vi = 1/(n*p_raw*(1-p_raw)); }
                        dv=e/n;
                    } else if(type==='survival'){
                        const hr=Engine.toNum(r.hr), l=Engine.toNum(r.lci), u=Engine.toNum(r.uci);
                        if(hr<=0)return null; es=Math.log(hr); const se=(Math.log(u)-Math.log(l))/3.92; vi=se*se; dv=hr;
                    }
                    if(!Engine.validNum(es)||!Engine.validNum(vi))return null;
                    const se=Math.sqrt(vi);
                    if(Metric.isLogRatio(met)||met==='Prop'){
                        if(met==='Prop'){ lo=Math.exp(es-1.96*se)/(1+Math.exp(es-1.96*se)); hi=Math.exp(es+1.96*se)/(1+Math.exp(es+1.96*se)); }
                        else { lo=Math.exp(es-1.96*se); hi=Math.exp(es+1.96*se); }
                    } else { lo=es-1.96*se; hi=es+1.96*se; }
                    return {id:r.id||'Study', year:Engine.toNum(r.year), group:r.group, es, vi, se, displayVal:dv, lo, hi, raw:r, metric:met, rowIdx:i, excluded:!!r.excluded};
                }).filter(x=>x);
            },
            pool:(eff,mod,meth,est,hk)=>{
                const act=eff.filter(e=>!e.excluded); if(!act.length)return null;
                let sw=0,swy=0; act.forEach(e=>{const w=1/e.vi;sw+=w;swy+=w*e.es;});
                const mu=swy/sw, Q=act.reduce((a,b)=>a+(1/b.vi)*(b.es-mu)**2,0), df=Math.max(0,act.length-1);
                const I2=Math.max(0,(Q-df)/Q)*100, t2=Engine.solveTau2(act,est,mod);
                let pes=0,pse=0;
                if(meth==='mh' && ['RR','OR'].includes(act[0].metric)){
                    let num=0,den=0; act.forEach(e=>{ const n1=parseFloat(e.raw.n1),n2=parseFloat(e.raw.n2),ev1=parseFloat(e.raw.e1),ev2=parseFloat(e.raw.e2), N=n1+n2;
                    if(act[0].metric==='RR'){num+=(ev1*n2)/N;den+=(ev2*n1)/N;} else {num+=(ev1*(n2-ev2))/N;den+=(ev2*(n1-ev1))/N;} });
                    pes=Math.log(num/den); pse=Math.sqrt(1/sw);
                } else {
                    let s=0,sy=0; act.forEach(e=>{let w=1/(e.vi+t2);s+=w;sy+=w*e.es;}); pes=sy/s; pse=Math.sqrt(1/s);
                }
                if(hk && mod==='random' && act.length>1){
                    let s=0; act.forEach(e=>{let w=1/(e.vi+t2);s+=w;});
                    let vhk = act.reduce((a,e)=>a+(1/(e.vi+t2))*(e.es-pes)**2,0)/((act.length-1)*s);
                    pse=Math.sqrt(Math.max(pse*pse, vhk));
                }
                const z=hk?Engine.getT(df):1.96;
                let disp,lo,hi;
                const m=act[0].metric;
                if(Metric.isLogRatio(m)||m==='Prop'){
                    if(m==='Prop'){ disp=Math.exp(pes)/(1+Math.exp(pes)); lo=Math.exp(pes-z*pse)/(1+Math.exp(pes-z*pse)); hi=Math.exp(pes+z*pse)/(1+Math.exp(pes+z*pse)); }
                    else { disp=Math.exp(pes); lo=Math.exp(pes-z*pse); hi=Math.exp(pes+z*pse); }
                } else { disp=pes; lo=pes-z*pse; hi=pes+z*pse; }

                let pil=null,pih=null;
                if(mod==='random' && act.length>2){
                    const tp=Engine.getT(act.length-2), pw=tp*Math.sqrt(t2+pse*pse);
                    if(Metric.isLogRatio(m)||m==='Prop'){
                        if(m==='Prop'){ pil=Math.exp(pes-pw)/(1+Math.exp(pes-pw)); pih=Math.exp(pes+pw)/(1+Math.exp(pes+pw)); }
                        else { pil=Math.exp(pes-pw); pih=Math.exp(pes+pw); }
                    } else { pil=pes-pw; pih=pes+pw; }
                }
                const pv=Engine.getP(pes/pse);
                eff.forEach(e=>{
                    if (e.excluded) { e.w=0; e.imp=0; e.q=0; }
                    else {
                        e.w=1/(e.vi+t2);
                        e.imp=(e.es-pes)**2/(e.vi+t2);
                        e.q=(1/e.vi)*(e.es-mu)**2; // Baujat Q
                    }
                });
                return {es:pes, se:pse, I2, Q, tau2:t2, studies:eff, displayVal:disp, lo, hi, pi_lo:pil, pi_hi:pih, pVal:pv};
            },
            eggersTest: (effects) => {
                const active = effects.filter(e=>!e.excluded);
                if(active.length<3) return {p:null, intercept:0};
                const x=active.map(e=>1/e.se), y=active.map(e=>e.es/e.se);
                const n=x.length, sx=x.reduce((a,b)=>a+b,0), sy=y.reduce((a,b)=>a+b,0);
                const sxy=x.reduce((a,b,i)=>a+b*y[i],0), sxx=x.reduce((a,b)=>a+b*b,0);
                const slope = (n*sxy - sx*sy)/(n*sxx - sx*sx);
                const int = (sy - slope*sx)/n;
                const resid = y.map((yi,i)=>yi-(int+slope*x[i]));
                const sse = resid.reduce((a,b)=>a+b*b,0);
                const t = int / Math.sqrt((sse/(n-2))*(sxx/(n*sxx-sx*sx)));
                return { p: Engine.getP(t), intercept: int };
            }
        };

if (typeof module!=='undefined'&&module.exports){ module.exports = { Metric, Engine }; }
