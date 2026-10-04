// Render the phase paper with cached KaTeX, bundled marked/Playwright and Edge.
// Run with the existing Node runtime. Preview files stay in the OS temp folder.
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
const paper = path.resolve(here, '../可组合认知结构与有限维量子理论_阶段论文.md');
const source = fs.readFileSync(paper, 'utf8');
const sha = x => crypto.createHash('sha256').update(x).digest('hex');
const display = [...source.matchAll(/^\$\$\s*\r?\n([\s\S]*?)^\$\$\s*$/gm)].map(m=>m[1]);
if ((source.match(/\$\$/g)||[]).length !== display.length*2) throw Error('Unpaired display delimiters');
const options = {throwOnError:true, strict:'ignore', trust:false};
let index = 0;
let transformed = source.replace(/^\$\$\s*\r?\n([\s\S]*?)^\$\$\s*$/gm, (_,tex)=>
  '\n<div class="formula" data-index="'+(++index)+'">'+katex.renderToString(tex,{...options,displayMode:true})+'</div>\n');
const prose = source.replace(/^```[\s\S]*?^```\s*$/gm,'').replace(/^\$\$\s*\r?\n[\s\S]*?^\$\$\s*$/gm,'');
const inline = [...prose.matchAll(/(?<!\\)\$([^$\n]+)\$/g)].map(m=>m[1]);
if ((prose.match(/(?<!\\)\$/g)||[]).length !== inline.length*2) throw Error('Unpaired inline delimiters');
inline.forEach(tex=>katex.renderToString(tex,{...options,displayMode:false}));
transformed = transformed.replace(/(?<!\\)\$([^$\n]+)\$/g,(_,tex)=>katex.renderToString(tex,{...options,displayMode:false}));
const tags = [...source.matchAll(/\\tag\{(\d+)\}/g)].map(m=>Number(m[1]));
if (tags.some((n,i)=>n!==i+1)) throw Error('Equation tags are not sequential');
const css = fs.readFileSync(path.join(cache,'katex.min.css'),'utf8');
const html = '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'+
  '<base href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/"><style>'+css+'</style><style>'+ 
  'body{margin:0;background:#eef1f5;color:#162237;font:17px/1.85 "Microsoft YaHei",sans-serif}'+
  'main{box-sizing:border-box;max-width:1032px;margin:24px auto;padding:52px 56px;background:white}'+
  'h1{font-size:30px;line-height:1.4}h2{font-size:23px;margin-top:2em}h3{font-size:19px;margin-top:1.6em}'+
  'p{margin:1em 0}a{color:#245b8d;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.7}'+
  'td,th{border:1px solid #c8d1dd;padding:8px 10px;text-align:left}th{background:#eff3f8}'+
  'pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:16px;background:#f3f5f8;font-size:13px}'+
  'code{overflow-wrap:anywhere}li{margin:8px 0}.formula{margin:24px 0}.katex{font-size:1.06em}'+
  '</style></head><body><main>'+marked.parse(transformed)+'</main></body></html>';
const outdir = path.join(os.tmpdir(),'codex-cognition-stage2-paper');
fs.mkdirSync(outdir,{recursive:true});
fs.writeFileSync(path.join(outdir,'paper.html'),html);
const browser = await chromium.launch({channel:'msedge',headless:true});
let report;
try {
  const page = await browser.newPage({viewport:{width:1280,height:1100},deviceScaleFactor:1});
  await page.setContent(html,{waitUntil:'networkidle',timeout:30000});
  await page.evaluate(()=>document.fonts.ready);
  const layout = await page.evaluate(()=>({
    display_count:document.querySelectorAll('.katex-display').length,
    errors:document.querySelectorAll('.katex-error').length,
    page_overflow:document.documentElement.scrollWidth>innerWidth,
    formula_overflow:[...document.querySelectorAll('.formula')].map(e=>({formula:Number(e.dataset.index),overflow:e.scrollWidth>e.clientWidth+1})).filter(x=>x.overflow),
    table_overflow:[...document.querySelectorAll('table')].map((e,i)=>({table:i+1,overflow:e.getBoundingClientRect().right>document.querySelector('main').getBoundingClientRect().right-30})).filter(x=>x.overflow),
    font_errors:[...document.fonts].filter(f=>f.status==='error').map(f=>f.family),
    document_height:document.documentElement.scrollHeight,
  }));
  if(layout.display_count!==display.length||layout.errors||layout.page_overflow||layout.formula_overflow.length||layout.table_overflow.length||layout.font_errors.length) throw Error(JSON.stringify(layout));
  const shots = [];
  for(const [name, heading] of [['opening',null],['theorem','3 综合重建定理'],['instruments','7 全部测量、通道和有限仪器'],['dynamics','8 连续可逆演化的生成元']]) {
    await page.evaluate(text=>{const e=text?[...document.querySelectorAll('h2')].find(e=>e.textContent===text):null;window.scrollTo(0,e?e.getBoundingClientRect().top+scrollY-20:0);},heading);
    const target=path.join(outdir,name+'.png');
    await page.screenshot({path:target});
    shots.push(target);
  }
  report={date:'2026-09-20',paper:path.basename(paper),paper_sha256:sha(fs.readFileSync(paper)),
    renderer:'KaTeX '+katex.version,renderer_sha256:sha(fs.readFileSync(path.join(cache,'katex.min.js'))),
    browser:'Headless Microsoft Edge via bundled Playwright',display_formulas:display.length,inline_formulas:inline.length,
    equation_tags:tags,render_layout:layout,
    scope:'All formulas parsed; complete paper rendered and programmatically checked for overflow and font errors. Representative screenshots saved for visual inspection.',
    temporary_preview_directory:outdir,representative_screenshots:shots,all_passed:true};
} finally {await browser.close();}
fs.writeFileSync(path.join(here,'stage2_math_review_checks.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
