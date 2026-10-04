// Parse and render the two closing notes and paper addendum without changing them.
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import crypto from 'node:crypto';
import {createRequire} from 'node:module';
import {fileURLToPath, pathToFileURL} from 'node:url';

const require = createRequire(import.meta.url);
const here = path.dirname(fileURLToPath(import.meta.url));
const runtime = path.join(os.homedir(), '.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules');
const cache = path.join(os.tmpdir(), 'codex-cognition-math-review');
const katex = require(path.join(cache, 'katex.min.js'));
const {marked} = await import(pathToFileURL(path.join(runtime, 'marked/lib/marked.esm.js')));
const {chromium} = require(path.join(runtime, 'playwright'));
const sha = x => crypto.createHash('sha256').update(x).digest('hex');
const displayPattern = /^\$\$\s*\r?\n([\s\S]*?)^\$\$\s*$/gm;
const options = {throwOnError:true, strict:'ignore', trust:false};
function render(tex, displayMode) {
  try { return katex.renderToString(tex,{...options,displayMode}); }
  catch (error) { throw new Error(error.message); }
}
const css = fs.readFileSync(path.join(cache,'katex.min.css'),'utf8');
const outdir = path.join(os.tmpdir(),'codex-cognition-stage2-completion');
fs.mkdirSync(outdir,{recursive:true});
const browser = await chromium.launch({channel:'msedge',headless:true});
const records = [];
try {
  for (const name of ['research_note_229.md','research_note_230.md','STAGE2_ADDENDUM.md']) {
    const source = fs.readFileSync(path.join(here,name),'utf8');
    const display = [...source.matchAll(displayPattern)].map(m=>m[1]);
    if ((source.match(/\$\$/g)||[]).length !== display.length*2) throw Error(name+': unpaired display delimiters');
    let index = 0;
    const transformed = source.replace(displayPattern,(_,tex)=>
      '\n<div class="formula" data-index="'+(++index)+'">'+render(tex,true)+'</div>\n');
    const prose = source.replace(displayPattern,'').replace(/^```[\s\S]*?^```\s*$/gm,'');
    const inline = [...prose.matchAll(/(?<!\\)\$([^$\n]+)\$/g)].map(m=>m[1]);
    if ((prose.match(/(?<!\\)\$/g)||[]).length!==inline.length*2) throw Error(name+': unpaired inline delimiters');
    inline.forEach(tex=>render(tex,false));
    const tags = [...source.matchAll(/\\tag\{(\d+)\}/g)].map(m=>Number(m[1]));
    if(tags.some((n,i)=>n!==i+1)) throw Error(name+': nonsequential equation tags');
    const html = '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'+
      '<base href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/"><style>'+css+'</style><style>'+
      'body{margin:0;background:#edf1f6;color:#19263b;font:17px/1.85 "Microsoft YaHei",sans-serif}'+
      'main{box-sizing:border-box;max-width:1080px;margin:24px auto;padding:44px 54px;background:white}'+
      'h1{font-size:29px;line-height:1.5}h2{font-size:23px;margin-top:1.8em}h3{font-size:19px}'+
      'a{color:#245b8d;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px}'+
      'td,th{border:1px solid #c8d1dd;padding:7px 9px;text-align:left}th{background:#eff3f8}'+
      'pre{white-space:pre-wrap;overflow-wrap:anywhere}li{margin:7px 0}.formula{margin:24px 0}.katex{font-size:1.02em}'+
      '</style></head><body><main>'+marked.parse(transformed)+'</main></body></html>';
    fs.writeFileSync(path.join(outdir,name+'.html'),html);
    const page = await browser.newPage({viewport:{width:1280,height:1100}});
    await page.setContent(html,{waitUntil:'networkidle',timeout:30000});
    await page.evaluate(()=>document.fonts.ready);
    const layout = await page.evaluate(()=>({
      display_count:document.querySelectorAll('.katex-display').length,
      errors:document.querySelectorAll('.katex-error').length,
      page_overflow:document.documentElement.scrollWidth>innerWidth,
      formula_overflow:[...document.querySelectorAll('.formula')].filter(e=>e.scrollWidth>e.clientWidth+1).map(e=>Number(e.dataset.index)),
      table_overflow:[...document.querySelectorAll('table')].some(e=>e.getBoundingClientRect().right>document.querySelector('main').getBoundingClientRect().right-30),
      font_errors:[...document.fonts].filter(f=>f.status==='error').map(f=>f.family),
    }));
    if(layout.display_count!==display.length||layout.errors||layout.page_overflow||layout.formula_overflow.length||layout.table_overflow||layout.font_errors.length) throw Error(name+': '+JSON.stringify(layout));
    const screenshot = path.join(outdir,name+'.png');
    await page.screenshot({path:screenshot});
    records.push({path:name,sha256:sha(fs.readFileSync(path.join(here,name))),display_formulas:display.length,inline_formulas:inline.length,equation_tags:tags,layout,screenshot});
    await page.close();
  }
} finally { await browser.close(); }
const report = {date:'2026-09-21',renderer:'KaTeX '+katex.version,browser:'Headless Edge via bundled Playwright',
  scope:'All three documents parsed and rendered; every formula and table checked for overflow. Opening screenshots saved for visual inspection.',
  documents:records,all_passed:true};
fs.writeFileSync(path.join(here,'completion_math_checks.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
