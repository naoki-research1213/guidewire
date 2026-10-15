"""JSON設定の読込みと、SI長さから既存計算座標への変換だけを担当する。"""
from pathlib import Path
import json, math

def load_parameters(path):
    p=json.loads(Path(path).read_text(encoding='utf-8'))
    if p['schema_version'] != 1: raise ValueError('Unsupported parameter schema')
    scale=float(p['internal_length_m'])
    if not math.isfinite(scale) or scale<=0: raise ValueError('internal_length_m must be positive')
    g=p['geometry_m']; m=p['mechanics_code']; n=p['numerics']; drive=p['input']; contact=p['contact']
    if 'material_si' in p:
        material=p['material_si']; mass=float(p['internal_mass_kg'])
        if not math.isfinite(mass) or mass<=0: raise ValueError('internal_mass_kg must be positive')
        # 時間単位は1秒。単位換算であり、原力学処理の置換ではない。
        m['young_modulus']=float(material['young_modulus_pa'])*scale/mass
        m['mass_density']=float(material['density_kg_m3'])*scale**3/mass
        m['poisson_ratio']=float(material['poisson_ratio_assumed'])
        m['gravity_y_code']=float(material['gravity_y_m_s2'])/scale
        if 'boundary_si' in p:
            m['support_stiffness']=float(p['boundary_si']['translation_stiffness_n_m'])/mass
            m['support_angular_stiffness']=float(p['boundary_si']['rotation_stiffness_n_m_rad'])/(mass*scale**2)
        p['force_scale_n_per_code']=mass*scale
    values=[n['dt_s'],g['radius']/scale,m['young_modulus'],m['mass_density'],m['poisson_ratio'],
        g['straight_length']/scale,g['tip_length']/scale,g['spire_diameter']/scale,
        drive['insertion_speed_m_s']/scale,drive['rotation_speed_rad_s'],m['rayleigh_mass'],m['rayleigh_stiffness'],
        m['support_stiffness'],m['support_angular_stiffness'],contact['alarm_distance_m']/scale,
        contact['distance_m']/scale,n['contact_tolerance_code'],n['contact_max_iterations'],
        m['gravity_y_code'],g['redistribution_threshold']/scale,contact['angle_cone'],contact['cone_factor'],contact['friction']]
    values=list(map(float,values))
    if not all(math.isfinite(x) for x in values): raise ValueError('Non-finite scene parameter')
    if values[0]<=0 or values[1]<=0 or values[2]<=0 or values[3]<=0: raise ValueError('Positive dt/radius/E/density required')
    if not 0<=values[4]<.5: raise ValueError('poisson_ratio must be in [0, 0.5)')
    if contact['distance_m']<0 or contact['alarm_distance_m']<=contact['distance_m']: raise ValueError('alarm distance must exceed contact distance')
    return p,values
