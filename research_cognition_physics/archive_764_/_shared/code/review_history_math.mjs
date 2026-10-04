// Render the new note selected by argument 241 or 242; old reports stay unchanged.
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import crypto from 'node:crypto';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const require=createRequire(import.meta.url);
const here=path.dirname(fileURLToPath(import.meta.url));
const round=Number(process.argv[2]);
if(![241,242].includes(round)) throw Error('Specify 241 or 242');
const round232=false;
const round231=true;
const base=here;
const names=['research_note_'+round+'.md'];
const runtime=path.join(os.homedir(),'.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules');
const cache=path.join(os.tmpdir(),'codex-cognition-math-review');
const katex=require(path.join(cache,'katex.min.js'));
const {marked}=await import(pathToFileURL(path.join(runtime,'marked/lib/marked.esm.js')));
const {chromium}=require(path.join(runtime,'playwright'));
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const displayPattern=/^\$\$\s*\r?\n([\s\S]*?)^\$\$\s*$/gm;
const inlinePattern=/(?<!\\)\$([^$\n]+)\$/g;
function render(tex,displayMode) {
  try { return katex.renderToString(tex,{displayMode,throwOnError:true,strict:'ignore',trust:false}); }
  catch(e) { throw Error(e.message); }
}
const css=fs.readFileSync(path.join(cache,'katex.min.css'),'utf8');
const outdir=path.join(os.tmpdir(),'codex-cognition-round'+round);
fs.mkdirSync(outdir,{recursive:true});
const browser=await chromium.launch({channel:'msedge',headless:true});
const records=[];
try {
  for (const name of names) {
    const filename=path.resolve(base,name);
    const source=fs.readFileSync(filename,'utf8');
    const display=[...source.matchAll(displayPattern)];
    if((source.match(/\$\$/g)||[]).length!==2*display.length) throw Error(name+': unpaired display delimiters');
    const prose=source.replace(displayPattern,'').replace(/^```[\s\S]*?^```\s*$/gm,'');
    const inline=[...prose.matchAll(inlinePattern)];
    if((prose.match(/(?<!\\)\$/g)||[]).length!==2*inline.length) throw Error(name+': unpaired inline delimiters');
    inline.forEach(m=>render(m[1],false));
    const tags=[...source.matchAll(/\\tag\{(\d+)\}/g)].map(m=>Number(m[1]));
    if(tags.some((v,i)=>v!==i+1)) throw Error(name+': nonsequential tags');
    // Replace formulas by tokens before Markdown parsing, so KaTeX HTML is not
    // reinterpreted as inline Markdown or another formula.
    const fragments=[];
    const token=html=>{fragments.push(html);return 'MATHTOKEN'+(fragments.length-1)+'ENDTOKEN';};
    let index=0;
    let transformed=source.replace(displayPattern,(_,tex)=>
      '\n\n'+token('<div class="formula" data-index="'+(++index)+'">'+render(tex,true)+'</div>')+'\n\n');
    transformed=transformed.replace(inlinePattern,(_,tex)=>token(render(tex,false)));
    let body=marked.parse(transformed)
      .replace(/<p>(MATHTOKEN\d+ENDTOKEN)<\/p>/g,(_,t)=>t)
      .replace(/MATHTOKEN(\d+)ENDTOKEN/g,(_,i)=>fragments[Number(i)]);
    const html='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'+
      '<base href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/"><style>'+css+'</style><style>'+
      'body{margin:0;background:#eef1f5;color:#162237;font:17px/1.85 "Microsoft YaHei",sans-serif}'+
      'main{box-sizing:border-box;max-width:1080px;margin:24px auto;padding:48px 54px;background:white}'+
      'h1{font-size:29px;line-height:1.5}h2{font-size:23px;margin-top:1.8em}h3{font-size:19px}'+
      'a{color:#245b8d;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.7}'+
      'td,th{border:1px solid #c8d1dd;padding:8px 10px;text-align:left}th{background:#eff3f8}'+
      'pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:14px;background:#f3f5f8;font-size:13px}'+
      'code{overflow-wrap:anywhere}li{margin:8px 0}.formula{margin:24px 0}.katex{font-size:1.04em}'+
      '</style></head><body><main>'+body+'</main></body></html>';
    const stem=path.basename(name,'.md');
    fs.writeFileSync(path.join(outdir,stem+'.html'),html);
    const page=await browser.newPage({viewport:{width:1280,height:1100}});
    await page.setContent(html,{waitUntil:'networkidle',timeout:30000});
    await page.evaluate(()=>document.fonts.ready);
    const layout=await page.evaluate(()=>({
      display_count:document.querySelectorAll('.katex-display').length,
      inline_count:document.querySelectorAll('.katex').length-document.querySelectorAll('.katex-display').length,
      errors:document.querySelectorAll('.katex-error').length,
      page_overflow:document.documentElement.scrollWidth>innerWidth,
      formula_overflow:[...document.querySelectorAll('.formula')].filter(e=>e.scrollWidth>e.clientWidth+1).map(e=>Number(e.dataset.index)),
      table_overflow:[...document.querySelectorAll('table')].some(e=>e.getBoundingClientRect().right>document.querySelector('main').getBoundingClientRect().right-30),
      font_errors:[...document.fonts].filter(f=>f.status==='error').map(f=>f.family)
    }));
    if(layout.display_count!==display.length||layout.inline_count!==inline.length||layout.errors||
       layout.page_overflow||layout.formula_overflow.length||layout.table_overflow||layout.font_errors.length)
      throw Error(name+': '+JSON.stringify(layout));
    const headings=round===241?[null,'3. 四条','4. 四种']:[null,'2. 相同','3. 有限'];
    const screenshots=[];
    for(const [i,prefix] of headings.entries()) {
      const found=await page.evaluate(prefix=>{
        const e=prefix?[...document.querySelectorAll('h2,h3')].find(e=>e.textContent.startsWith(prefix)):null;
        if(prefix&&!e) return false;
        scrollTo(0,e?e.getBoundingClientRect().top+scrollY-20:0);
        return true;
      },prefix);
      if(!found) continue;
      await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
      const target=path.join(outdir,stem+'-'+i+'.png');
      await page.screenshot({path:target});screenshots.push(target);
    }
    records.push({path:name,sha256:sha(fs.readFileSync(filename)),display_formulas:display.length,
      inline_formulas:inline.length,equation_tags:tags,layout,screenshots});
    await page.close();
  }
} finally {await browser.close();}
const report={date:'2026-09-22',renderer:'KaTeX '+katex.version,
  renderer_sha256:sha(fs.readFileSync(path.join(cache,'katex.min.js'))),
  browser:'Headless Edge via bundled Playwright',
  scope:'All display and inline formulas parsed and rendered; formula/table overflow and fonts checked. Screenshots for representative visual inspection.',
  documents:records,all_passed:true};
fs.writeFileSync(path.join(base,'round'+round+'_math_checks.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
