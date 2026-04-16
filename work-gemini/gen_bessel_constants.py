import re
from mpmath import mp, mpf, pi

mp.dps = 60

def to_dd(v):
    if isinstance(v, str):
        v = mpf(v)
    hi = float(v)
    lo = float(v - mpf(hi))
    return hi, lo

def extract(path):
    with open(path, 'r') as f:
        content = f.read()
    pattern = re.compile(r'static const __float128\s+(\w+)\[[^\]]*\]\s*=\s*\{([^\}]+)\};', re.MULTILINE)
    res = {}
    for match in pattern.finditer(content):
        name = match.group(1)
        vals_str = match.group(2)
        vals = []
        for v in vals_str.split(','):
            v = re.sub(r'/\*.*?\*/', '', v, flags=re.DOTALL)
            v = v.strip().replace('Q', '').replace('L', '')
            if v:
                vals.append(v)
        res[name] = vals
    return res

j0_data = extract('external/libquadmath/math/j0q.c')
j1_data = extract('external/libquadmath/math/j1q.c')

# Trig and other scalars
constants = {
    'INV_SQRT2': '0.7071067811865475244008443621048490392848',
    'TWO_OVER_PI': '0.6366197723675813430755350534900574481378',
    'U0': '-7.3804295108687225274343927948483016310862e-02',
    'PI_QUARTER': pi / 4,
    'THREE_PI_QUARTER': 3 * pi / 4
}

with open('work-gemini/bessel_constants.hh', 'w') as f:
    f.write('#pragma once\n\nnamespace bessel_improved {\n\n')
    for name, val in constants.items():
        hi, lo = to_dd(val)
        f.write(f'static const double {name}_hi = {repr(hi)};\n')
        f.write(f'static const double {name}_lo = {repr(lo)};\n\n')
    
    for prefix, data in [('j0', j0_data), ('j1', j1_data)]:
        for name, vals in data.items():
            full_name = f'{prefix}_{name}'
            if prefix == 'j1':
                if 'J0' in name: full_name = full_name.replace('J0', 'J1')
                if 'Y0' in name: full_name = full_name.replace('Y0', 'Y1')
            
            f.write(f'static const double {full_name}_hi[] = {{\n')
            for v in vals:
                hi, _ = to_dd(v)
                f.write(f'  {repr(hi)},\n')
            f.write('};\n')
            f.write(f'static const double {full_name}_lo[] = {{\n')
            for v in vals:
                _, lo = to_dd(v)
                f.write(f'  {repr(lo)},\n')
            f.write('};\n\n')
    f.write('} // namespace bessel_improved\n')
