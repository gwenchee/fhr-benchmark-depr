""" Constants and functions for FHR benchmark geometry generation

This script contains the values for constants and functions used to 
set up the FHR benchmark geometry. 

"""

from numpy import sin, cos, tan, pi
import openmc

###############################################################################
#                               Constants 
###############################################################################

H_side = 22.5/sin(pi/3)
P_len = 23.1 # plank length
P_D_jut = 2-1.4948
P_D_jut_hyp = P_D_jut/sin(pi/3)
P_D_jut_adj = P_D_jut/sin(pi/3)
P_small_gap = 0.35
P_A1_height = 2.55
P_A1_adj = P_A1_height/tan(pi/3)
P_big_gap = 0.7
P_A2_hyp =  P_A1_height/sin(pi/3)
P_big_gap_A2_hyp = P_big_gap/sin(pi/3)
P_A3_hyp = P_A1_height/sin(pi/3)
P_big_gap_A3_hyp = P_big_gap/sin(pi/3)
D_to_center = 2
D_to_center_width = D_to_center*tan(pi/6)
D_A1_width = P_len - 2*(P_D_jut)
D_A1_height = 19.5
D_A1_adj = D_A1_height/tan(pi/3)
T_pitch = 0.09266
T_r1 = 2135e-5
T_r2 = 3135e-5
T_r3 = 3485e-5
T_r4 = 3835e-5
T_r5 = 4235e-5
F_protect_gap = 0.1
F_width = T_pitch*(4)
F_len = T_pitch*(210)
F_A1_D_gap = (D_A1_width-F_len)/2
F_F_gap = P_A1_height - 2*F_width -2*F_protect_gap
F_F_gap_adj = F_F_gap/tan(pi/3)
F_A1_width_adj = F_width/tan(pi/3)
F_F_gap_A2_hyp = F_F_gap/sin(pi/3)
F_A2_width_hyp = F_width/sin(pi/3)
F_F_gap_A3_hyp = F_F_gap/sin(pi/3)
F_F_gap_A3_adj = F_F_gap_A3_hyp*cos(pi/3)
F_F_gap_A3_opp = F_F_gap_A3_hyp*sin(pi/3)
F_A3_width_adj = F_width*cos(pi/3)
F_A3_width_opp = F_width*sin(pi/3)
S_S_gap = 14
S_A1_D_gap = (D_A1_width-S_S_gap)/2
S_large_r = 0.7
S_small_r = 0.35
CS_l = 10.38
CS_w = 1.76
CA_l = 10
CA_w = 1 
z_thickness = 101 # no. of triso particle thickness, must be odd 

###############################################################################
#                                  Functions
###############################################################################

def plane(m,x,y): 
    """ Creates openmc plane

    Parameters
    ----------
    m: float 
        gradient of plane 
    x: float 
        x-intercept 
    y: float
        y-intercept 

    Returns
    -------
    openmc.Plane 
    """

    return openmc.Plane(a=-m,b=1,d=-m*x+y)

def region_maker(area,area_type):
    """ Creates 2D openmc.region.Intersection on x-y plane 

    Parameters
    ----------
    area: str 
        This refers to the 3 diamond sections in the FHR geometry. Options are:
        'A1', 'A2', 'A3'
    area_type: str
        This refers to the various areas within each quadrant. 
        Diamond Plank Area: 'D', Plank Area: 'P', Fuel Area: 'F' 
        Spacers: 'S', Control Rod Slot: 'CS', Control Rod Arm: 'CA' 

    Returns
    -------
    openmc.region.Intersection
    """

    if area in ['A1','A3']:
        if V[area][area_type]['L']['m'] == 0.0 and V[area][area_type]['R']['m'] == 0.0: 
            region = -plane(V[area][area_type]['T']['m'],V[area][area_type]['T']['x'],V[area][area_type]['T']['y']) &\
                     +plane(V[area][area_type]['B']['m'],V[area][area_type]['B']['x'],V[area][area_type]['B']['y']) &\
                     +openmc.XPlane(x0=V[area][area_type]['L']['x']) &\
                     -openmc.XPlane(x0=V[area][area_type]['R']['x'])     
        else: 
            region = -plane(V[area][area_type]['T']['m'],V[area][area_type]['T']['x'],V[area][area_type]['T']['y']) &\
                     +plane(V[area][area_type]['B']['m'],V[area][area_type]['B']['x'],V[area][area_type]['B']['y']) &\
                     +plane(V[area][area_type]['L']['m'],V[area][area_type]['L']['x'],V[area][area_type]['L']['y']) &\
                     -plane(V[area][area_type]['R']['m'],V[area][area_type]['R']['x'],V[area][area_type]['R']['y']) 
        
    elif area in ['A2']:
        region = -plane(V[area][area_type]['T']['m'],V[area][area_type]['T']['x'],V[area][area_type]['T']['y']) &\
                 +plane(V[area][area_type]['B']['m'],V[area][area_type]['B']['x'],V[area][area_type]['B']['y']) &\
                 -plane(V[area][area_type]['L']['m'],V[area][area_type]['L']['x'],V[area][area_type]['L']['y']) &\
                 +plane(V[area][area_type]['R']['m'],V[area][area_type]['R']['x'],V[area][area_type]['R']['y']) 
    else: 
        raise Exception('Your region type has yet to be defined.')

    return region 

def rx(x_i,y_i,t):
    """ Rotates x on x-y plane counterclockwise 
    through a specified angle with respect to the x axis about the 
    origin of a 2D cartesian coordinate system. 

    Parameters
    ----------
    x_i float 
        x value to rotate 
    y_i: float 
        y value to rotate
    t: float
        angle to rotate about origin 

    Returns
    -------
    rotated x value 
    """

    return x_i*cos(t)-y_i*sin(t)

def ry(x_i,y_i,t):
    """ Rotates y on x-y plane counterclockwise 
    through a specified angle with respect to the x axis about the 
    origin of a 2D cartesian coordinate system. 

    Parameters
    ----------
    x_i float 
        x value to rotate 
    y_i: float 
        y value to rotate
    t: float
        angle to rotate about origin 

    Returns
    -------
    rotated y value 
    """
    return x_i*sin(t)+y_i*cos(t)

###############################################################################
#                             Geometry for Planes 
# The V (values) dictionary organizes the x-intercept, y-intercept, and 
# gradient values required to create the planes that set up the regions for 
# various areas of the FHR geometry. 
# Example V format: V[level1][level2][level3][level4]
# level1: 3 diamond sections in the FHR geometry. Options are: 'A1', 'A2', 'A3'
# level2: various areas within each quadrant. Options are:
#         Diamond Plank Area: 'D', Plank Area: 'P', Fuel Area: 'F' 
#         Spacers: 'S', Control Rod Slot: 'CS', Control Rod Arm: 'CA' 
# level3: plane position relative to area. Options are: 
#         Top: 'T', Bottom: 'B', Left: 'L', Right: 'R'
# level4: gradient, x-intercept, or y-intercept value. Options are: 
#         gradient: 'm', x-intercept: 'x', y-intercept: 'y'
###############################################################################

m1 = -D_A1_height/D_A1_adj
m2 = D_A1_width*sin(pi/3)/(D_A1_width*cos(pi/3))

V = {'A1':{},'A2':{},'A3':{}}

V['A1'] = {'D':{'T':{},'B':{},'L':{},'R':{}},
           'P':{'T':{},'B':{},'L':{},'R':{}},
           'F':{'T':{},'B':{},'L':{},'R':{}},
           'S':{'C':{},'Cb':{}},
           'CS':{'T':{},'B':{},'L':{},'R':{}},
           'CA':{'T':{},'B':{},'L':{},'R':{}}}
V['A2'] = {'D':{'T':{},'B':{},'L':{},'R':{}},
           'P':{'T':{},'B':{},'L':{},'R':{}},
           'F':{'T':{},'B':{},'L':{},'R':{}},
           'S':{'C':{},'Cb':{}},
           'CS':{'T':{},'B':{},'L':{},'R':{}},
           'CA':{'T':{},'B':{},'L':{},'R':{}}}
V['A3'] = {'D':{'T':{},'B':{},'L':{},'R':{}},
           'P':{'T':{},'B':{},'L':{},'R':{}},
           'F':{'T':{},'B':{},'L':{},'R':{}},
           'S':{'C':{},'Cb':{}},
           'CS':{'T':{},'B':{},'L':{},'R':{}},
           'CA':{'T':{},'B':{},'L':{},'R':{}}}

V['A1']['D']['T'] = {'m':0, 'x':0, 'y':-D_to_center}
V['A1']['D']['B'] = {'m':0, 'x':0, 'y':-D_to_center-D_A1_height}
V['A1']['D']['R'] = {'m':m1, 'x':-D_to_center_width, 'y':-D_to_center}
V['A1']['D']['L'] = {'m':m1, 'x':V['A1']['D']['R']['x']-D_A1_width, 'y':-D_to_center}

V['A1']['P']['T'] = {'m':0, 'x':0, 'y':V['A1']['D']['T']['y']-P_small_gap}
V['A1']['P']['B'] = {'m':0, 'x':0, 'y':V['A1']['P']['T']['y']-P_A1_height}
V['A1']['P']['R'] = {'m':m1, 'x':V['A1']['D']['R']['x']+P_small_gap*tan(pi/6)+P_D_jut_hyp, 'y':V['A1']['D']['R']['y']-P_small_gap}
V['A1']['P']['L'] = {'m':m1, 'x':V['A1']['D']['L']['x']+P_small_gap*tan(pi/6)-P_D_jut_hyp, 'y':V['A1']['D']['L']['y']-P_small_gap}

V['A1']['F']['T'] = {'m':0, 'x':V['A1']['D']['R']['x']-F_A1_D_gap, 'y':V['A1']['P']['T']['y']-F_protect_gap}
V['A1']['F']['B'] = {'m':0, 'x':V['A1']['F']['T']['x'], 'y':V['A1']['F']['T']['y']-F_width}
V['A1']['F']['R'] = {'m':0, 'x':V['A1']['D']['R']['x']-F_A1_D_gap, 'y':V['A1']['F']['T']['y']}
V['A1']['F']['L'] = {'m':0, 'x':V['A1']['D']['L']['x']+F_A1_D_gap, 'y':V['A1']['F']['T']['y']}

V['A1']['S']['C'] = {'x0':-D_to_center_width-S_A1_D_gap,'y0':-D_to_center-P_small_gap}
V['A1']['S']['Cb'] = {'x0':-D_to_center_width-S_A1_D_gap+P_A1_adj,'y0':-D_to_center-P_small_gap-P_A1_height}

V['A1']['CS']['T'] = {'m':0, 'x':0, 'y':CS_w/2}
V['A1']['CS']['B'] = {'m':0, 'x':0, 'y':-CS_w/2}
V['A1']['CS']['R'] = {'m':0, 'x':0, 'y':0}
V['A1']['CS']['L'] = {'m':0, 'x':-CS_l, 'y':0}

V['A1']['CA']['T'] = {'m':0, 'x':0, 'y':CA_w/2}
V['A1']['CA']['B'] = {'m':0, 'x':0, 'y':-CA_w/2}
V['A1']['CA']['R'] = {'m':0, 'x':0, 'y':0}
V['A1']['CA']['L'] = {'m':0, 'x':-CA_l, 'y':0}

A2_t = -pi/3*2
V['A2']['D']['T'] = {'m':0, 'x':rx(V['A1']['D']['L']['x'],V['A1']['D']['L']['y'],A2_t), 'y':ry(V['A1']['D']['L']['x'],V['A1']['D']['L']['y'],A2_t)}
V['A2']['D']['B'] = {'m':0, 'x':rx(V['A1']['D']['R']['x'],V['A1']['D']['R']['y'],A2_t), 'y':ry(V['A1']['D']['R']['x'],V['A1']['D']['R']['y'],A2_t)}
V['A2']['D']['R'] = {'m':m2, 'x':rx(V['A1']['D']['T']['x'],V['A1']['D']['T']['y'],A2_t), 'y':ry(V['A1']['D']['T']['x'],V['A1']['D']['T']['y'],A2_t)}
V['A2']['D']['L'] = {'m':m2, 'x':rx(V['A1']['D']['B']['x'],V['A1']['D']['B']['y'],A2_t), 'y':ry(V['A1']['D']['B']['x'],V['A1']['D']['B']['y'],A2_t)}

V['A2']['P']['T'] = {'m':0, 'x':rx(V['A1']['P']['L']['x'],V['A1']['P']['L']['y'],A2_t), 'y':ry(V['A1']['P']['L']['x'],V['A1']['P']['L']['y'],A2_t)}
V['A2']['P']['B'] = {'m':0, 'x':rx(V['A1']['P']['R']['x'],V['A1']['P']['R']['y'],A2_t), 'y':ry(V['A1']['P']['R']['x'],V['A1']['P']['R']['y'],A2_t)}
V['A2']['P']['R'] = {'m':m2, 'x':rx(V['A1']['P']['T']['x'],V['A1']['P']['T']['y'],A2_t), 'y':ry(V['A1']['P']['T']['x'],V['A1']['P']['T']['y'],A2_t)}
V['A2']['P']['L'] = {'m':m2, 'x':rx(V['A1']['P']['B']['x'],V['A1']['P']['B']['y'],A2_t), 'y':ry(V['A1']['P']['B']['x'],V['A1']['P']['B']['y'],A2_t)}

V['A2']['F']['T'] = {'m':-1/m2, 'x':rx(V['A1']['F']['L']['x'],V['A1']['F']['L']['y'],A2_t), 'y':ry(V['A1']['F']['L']['x'],V['A1']['F']['L']['y'],A2_t)}
V['A2']['F']['B'] = {'m':-1/m2, 'x':rx(V['A1']['F']['R']['x'],V['A1']['F']['R']['y'],A2_t), 'y':ry(V['A1']['F']['R']['x'],V['A1']['F']['R']['y'],A2_t)}
V['A2']['F']['R'] = {'m':m2, 'x':rx(V['A1']['F']['T']['x'],V['A1']['F']['T']['y'],A2_t), 'y':ry(V['A1']['F']['T']['x'],V['A1']['F']['T']['y'],A2_t)}
V['A2']['F']['L'] = {'m':m2, 'x':rx(V['A1']['F']['B']['x'],V['A1']['F']['B']['y'],A2_t), 'y':ry(V['A1']['F']['B']['x'],V['A1']['F']['B']['y'],A2_t)}

V['A2']['S']['C'] = {'x0':rx(V['A1']['S']['C']['x0'],V['A1']['S']['C']['y0'],A2_t),'y0':ry(V['A1']['S']['C']['x0'],V['A1']['S']['C']['y0'],A2_t)}
V['A2']['S']['Cb'] = {'x0':rx(V['A1']['S']['Cb']['x0'],V['A1']['S']['Cb']['y0'],A2_t),'y0':ry(V['A1']['S']['Cb']['x0'],V['A1']['S']['Cb']['y0'],A2_t)}

V['A2']['CS']['T'] = {'m':-1/m2, 'x':rx(V['A1']['CS']['L']['x'],V['A1']['CS']['L']['y'],A2_t), 'y':ry(V['A1']['CS']['L']['x'],V['A1']['CS']['L']['y'],A2_t)}
V['A2']['CS']['B'] = {'m':-1/m2, 'x':rx(V['A1']['CS']['R']['x'],V['A1']['CS']['R']['y'],A2_t), 'y':ry(V['A1']['CS']['R']['x'],V['A1']['CS']['R']['y'],A2_t)}
V['A2']['CS']['R'] = {'m':m2, 'x':rx(V['A1']['CS']['T']['x'],V['A1']['CS']['T']['y'],A2_t), 'y':ry(V['A1']['CS']['T']['x'],V['A1']['CS']['T']['y'],A2_t)}
V['A2']['CS']['L'] = {'m':m2, 'x':rx(V['A1']['CS']['B']['x'],V['A1']['CS']['B']['y'],A2_t), 'y':ry(V['A1']['CS']['B']['x'],V['A1']['CS']['B']['y'],A2_t)}

V['A2']['CA']['T'] = {'m':-1/m2, 'x':rx(V['A1']['CA']['L']['x'],V['A1']['CA']['L']['y'],A2_t), 'y':ry(V['A1']['CA']['L']['x'],V['A1']['CA']['L']['y'],A2_t)}
V['A2']['CA']['B'] = {'m':-1/m2, 'x':rx(V['A1']['CA']['R']['x'],V['A1']['CA']['R']['y'],A2_t), 'y':ry(V['A1']['CA']['R']['x'],V['A1']['CA']['R']['y'],A2_t)}
V['A2']['CA']['R'] = {'m':m2, 'x':rx(V['A1']['CA']['T']['x'],V['A1']['CA']['T']['y'],A2_t), 'y':ry(V['A1']['CA']['T']['x'],V['A1']['CA']['T']['y'],A2_t)}
V['A2']['CA']['L'] = {'m':m2, 'x':rx(V['A1']['CA']['B']['x'],V['A1']['CA']['B']['y'],A2_t), 'y':ry(V['A1']['CA']['B']['x'],V['A1']['CA']['B']['y'],A2_t)}

A3_t = pi/3*2
V['A3']['D']['T'] = {'m':m2, 'x':rx(V['A1']['D']['R']['x'],V['A1']['D']['R']['y'],A3_t), 'y':ry(V['A1']['D']['R']['x'],V['A1']['D']['R']['y'],A3_t)}
V['A3']['D']['B'] = {'m':m2, 'x':rx(V['A1']['D']['L']['x'],V['A1']['D']['L']['y'],A3_t), 'y':ry(V['A1']['D']['L']['x'],V['A1']['D']['L']['y'],A3_t)}
V['A3']['D']['R'] = {'m':m1, 'x':rx(V['A1']['D']['B']['x'],V['A1']['D']['B']['y'],A3_t), 'y':ry(V['A1']['D']['B']['x'],V['A1']['D']['B']['y'],A3_t)}
V['A3']['D']['L'] = {'m':m1, 'x':rx(V['A1']['D']['T']['x'],V['A1']['D']['T']['y'],A3_t), 'y':ry(V['A1']['D']['T']['x'],V['A1']['D']['T']['y'],A3_t)}

V['A3']['P']['T'] = {'m':m2, 'x':rx(V['A1']['P']['R']['x'],V['A1']['P']['R']['y'],A3_t), 'y':ry(V['A1']['P']['R']['x'],V['A1']['P']['R']['y'],A3_t)}
V['A3']['P']['B'] = {'m':m2, 'x':rx(V['A1']['P']['L']['x'],V['A1']['P']['L']['y'],A3_t), 'y':ry(V['A1']['P']['L']['x'],V['A1']['P']['L']['y'],A3_t)}
V['A3']['P']['R'] = {'m':m1, 'x':rx(V['A1']['P']['B']['x'],V['A1']['P']['B']['y'],A3_t), 'y':ry(V['A1']['P']['B']['x'],V['A1']['P']['B']['y'],A3_t)}
V['A3']['P']['L'] = {'m':m1, 'x':rx(V['A1']['P']['T']['x'],V['A1']['P']['T']['y'],A3_t), 'y':ry(V['A1']['P']['T']['x'],V['A1']['P']['T']['y'],A3_t)}

V['A3']['F']['T'] = {'m':-1/m1, 'x':rx(V['A1']['F']['R']['x'],V['A1']['F']['R']['y'],A3_t), 'y':ry(V['A1']['F']['R']['x'],V['A1']['F']['R']['y'],A3_t)}
V['A3']['F']['B'] = {'m':-1/m1, 'x':rx(V['A1']['F']['L']['x'],V['A1']['F']['L']['y'],A3_t), 'y':ry(V['A1']['F']['L']['x'],V['A1']['F']['L']['y'],A3_t)}
V['A3']['F']['R'] = {'m':m1, 'x':rx(V['A1']['F']['B']['x'],V['A1']['F']['B']['y'],A3_t), 'y':ry(V['A1']['F']['B']['x'],V['A1']['F']['B']['y'],A3_t)}
V['A3']['F']['L'] = {'m':m1, 'x':rx(V['A1']['F']['T']['x'],V['A1']['F']['T']['y'],A3_t), 'y':ry(V['A1']['F']['T']['x'],V['A1']['F']['T']['y'],A3_t)}

V['A3']['S']['C'] = {'x0':rx(V['A1']['S']['C']['x0'],V['A1']['S']['C']['y0'],A3_t),'y0':ry(V['A1']['S']['C']['x0'],V['A1']['S']['C']['y0'],A3_t)}
V['A3']['S']['Cb'] = {'x0':rx(V['A1']['S']['Cb']['x0'],V['A1']['S']['Cb']['y0'],A3_t),'y0':ry(V['A1']['S']['Cb']['x0'],V['A1']['S']['Cb']['y0'],A3_t)}

V['A3']['CS']['T'] = {'m':-1/m1, 'x':rx(V['A1']['CS']['R']['x'],V['A1']['CS']['R']['y'],A3_t), 'y':ry(V['A1']['CA']['R']['x'],V['A1']['CA']['R']['y'],A3_t)}
V['A3']['CS']['B'] = {'m':-1/m1, 'x':rx(V['A1']['CS']['L']['x'],V['A1']['CS']['L']['y'],A3_t), 'y':ry(V['A1']['CA']['L']['x'],V['A1']['CA']['L']['y'],A3_t)}
V['A3']['CS']['R'] = {'m':m1, 'x':rx(V['A1']['CS']['B']['x'],V['A1']['CS']['B']['y'],A3_t), 'y':ry(V['A1']['CA']['B']['x'],V['A1']['CA']['B']['y'],A3_t)}
V['A3']['CS']['L'] = {'m':m1, 'x':rx(V['A1']['CS']['T']['x'],V['A1']['CS']['T']['y'],A3_t), 'y':ry(V['A1']['CA']['T']['x'],V['A1']['CA']['T']['y'],A3_t)}

V['A3']['CA']['T'] = {'m':-1/m1, 'x':rx(V['A1']['CA']['R']['x'],V['A1']['CA']['R']['y'],A3_t), 'y':ry(V['A1']['CA']['R']['x'],V['A1']['CA']['R']['y'],A3_t)}
V['A3']['CA']['B'] = {'m':-1/m1, 'x':rx(V['A1']['CA']['L']['x'],V['A1']['CA']['L']['y'],A3_t), 'y':ry(V['A1']['CA']['L']['x'],V['A1']['CA']['L']['y'],A3_t)}
V['A3']['CA']['R'] = {'m':m1, 'x':rx(V['A1']['CA']['B']['x'],V['A1']['CA']['B']['y'],A3_t), 'y':ry(V['A1']['CA']['B']['x'],V['A1']['CA']['B']['y'],A3_t)}
V['A3']['CA']['L'] = {'m':m1, 'x':rx(V['A1']['CA']['T']['x'],V['A1']['CA']['T']['y'],A3_t), 'y':ry(V['A1']['CA']['T']['x'],V['A1']['CA']['T']['y'],A3_t)}


###############################################################################
#                            Translation for Planes
###############################################################################
T = {'A1':{'P':{},'F':{},'S':{}},'A2':{'P':{},'F':{},'S':{}},'A3':{'P':{},'F':{},'S':{}}}

T['A1']['P'] = {'x':(P_big_gap+P_A1_height)/tan(pi/3),'y':-(P_big_gap+P_A1_height)}
T['A2']['P'] = {'x':-P_A2_hyp-P_big_gap_A2_hyp, 'y':0}
T['A3']['P'] = {'x':(P_A3_hyp+P_big_gap_A3_hyp)*cos(pi/3), 'y':(P_A3_hyp+P_big_gap_A3_hyp)*sin(pi/3)}

T['A1']['F'] = {'x':F_F_gap_adj+F_A1_width_adj, 'y':-F_F_gap-F_width}
T['A2']['F'] = {'x':-F_F_gap_A2_hyp-F_A2_width_hyp, 'y':0}
T['A3']['F'] = {'x':F_F_gap_A3_adj+F_A3_width_adj, 'y':F_F_gap_A3_opp+F_A3_width_opp}

T['A1']['S'] = {'x':-S_S_gap, 'y':0}
T['A2']['S'] = {'x':S_S_gap*cos(pi/3), 'y':S_S_gap*sin(pi/3)}
T['A3']['S'] = {'x':S_S_gap*cos(pi/3), 'y':-S_S_gap*sin(pi/3)}