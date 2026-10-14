"""Independent complex-step check; no import of author code."""
from pathlib import Path
from fractions import Fraction as F
import cmath,json,argparse,hashlib
HERE=Path(__file__).resolve().parent
def yu(z):
    f=1-z/2;b=1+z/2
    return z*cmath.sqrt(b)/f**1.5,1/(f**1.5*cmath.sqrt(b))
def run():
    x=-.5;h=1.e-20
    y,u=yu(x)
    yp=yu(x+1j*h)[0].imag/h
    up=yu(x+1j*h)[1].imag/h
    uy=up/yp
    uy_x=(1/(2*(1+(x+1j*h)/2))).imag/h
    uyy=uy_x/yp
    cs2=uy/(uy+2*y.real*uyy)
    assert abs(uy-2/3)<1.e-14 and abs(cs2-4/9)<1.e-14
    assert abs(yp-.33049457887636624)<1.e-14
    xq=F(-1,2)
    sound=(1+xq)/(1+xq*xq/2)
    assert sound==F(4,9)
    hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('README.md','check.py','results.json','review.md','independent_check.py')}
    return {'status':'PASS_WITH_SCOPE','Y_prime':yp,'U_Y':uy,'U_YY':uyy,
      'sound_speed_squared':cs2,'exact_sound_speed_squared':str(sound),
      'author_check_imported':False,'bound_files_sha256':hashes,
      'scope':'local classical algebra only; no original instrument or quantum process'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true');a=ap.parse_args()
    obj=run();p=HERE/'independent_results.json'
    if a.write_results:
        with p.open('x',encoding='utf-8') as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
    elif p.exists():assert obj==json.loads(p.read_text(encoding='utf-8'))
    print(json.dumps(obj,ensure_ascii=False,indent=2))
