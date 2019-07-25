#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy as np

import matplotlib.pyplot as plt


#===========================================================================================================
def print_text(text,window_title,plot_title,x=0,y=0,fontsize=12):
    """
    Print a text bloc on a figure
    """

    n = len(text)

    fig = plt.figure(figsize=(5,n/3))
    fig.canvas.set_window_title(window_title)
    fig.suptitle(plot_title, fontsize=14)

    n = len(text)

    plt.axis([0, 20, 0, n+2])
    plt.axis("off")

    for j in range(n):
        plt.text(x, n+1-(y+j), text[j], ha='left', fontsize=fontsize, rotation=0, wrap=True)

    plt.show()


#===========================================================================================================
def draw_3d_view(aircraft,window_title,plot_title):
    """
    Build a 3 viewx drawing of the airplane
    """

    fus_width = aircraft.fuselage.width
    fus_height = aircraft.fuselage.height
    fus_length = aircraft.fuselage.length

    htp_span = aircraft.horizontal_tail.span
    htp_dihedral = aircraft.horizontal_tail.dihedral
    htp_t_o_c = aircraft.horizontal_tail.t_o_c
    htp_x_axe = aircraft.horizontal_tail.x_axe
    htp_z_axe = aircraft.horizontal_tail.z_axe
    htp_c_axe = aircraft.horizontal_tail.c_axe
    htp_x_tip = aircraft.horizontal_tail.x_tip
    htp_z_tip = aircraft.horizontal_tail.z_tip
    htp_c_tip = aircraft.horizontal_tail.c_tip

    vtp_t_o_c = aircraft.vertical_tail.t_o_c
    vtp_x_root = aircraft.vertical_tail.x_root
    vtp_z_root = aircraft.vertical_tail.z_root
    vtp_c_root = aircraft.vertical_tail.c_root
    vtp_x_tip = aircraft.vertical_tail.x_tip
    vtp_z_tip = aircraft.vertical_tail.z_tip
    vtp_c_tip = aircraft.vertical_tail.c_tip

    wing_x_root = aircraft.wing.x_root
    wing_y_root = aircraft.wing.y_root
    wing_z_root = aircraft.wing.z_root
    wing_c_root = aircraft.wing.c_root
    wing_toc_r = aircraft.wing.t_o_c_r
    wing_x_kink = aircraft.wing.x_kink
    wing_y_kink = aircraft.wing.y_kink
    wing_z_kink = aircraft.wing.z_kink
    wing_c_kink = aircraft.wing.c_kink
    wing_toc_k = aircraft.wing.t_o_c_k
    wing_x_tip = aircraft.wing.x_tip
    wing_y_tip = aircraft.wing.y_tip
    wing_z_tip = aircraft.wing.z_tip
    wing_c_tip = aircraft.wing.c_tip
    wing_toc_t = aircraft.wing.t_o_c_t

    if (aircraft.propulsion.architecture=="TP"):
        nac_length = aircraft.turbofan_nacelle.length
        nac_height = aircraft.turbofan_nacelle.width
        nac_width = aircraft.turbofan_nacelle.width
        nac_x_ext = aircraft.turbofan_nacelle.x_ext
        nac_y_ext = aircraft.turbofan_nacelle.y_ext
        nac_z_ext = aircraft.turbofan_nacelle.z_ext
    elif (aircraft.propulsion.architecture=="TF"):
        nac_length = aircraft.turbofan_nacelle.length
        nac_height = aircraft.turbofan_nacelle.width
        nac_width = aircraft.turbofan_nacelle.width
        nac_x_ext = aircraft.turbofan_nacelle.x_ext
        nac_y_ext = aircraft.turbofan_nacelle.y_ext
        nac_z_ext = aircraft.turbofan_nacelle.z_ext
        if (aircraft.turbofan_engine.n_engine==4):
            nac_x_int = aircraft.turbofan_nacelle.x_int
            nac_y_int = aircraft.turbofan_nacelle.y_int
            nac_z_int = aircraft.turbofan_nacelle.z_int
    elif (aircraft.propulsion.architecture=="EF1"):
        nac_length = aircraft.electrofan_nacelle.length
        nac_height = aircraft.electrofan_nacelle.width
        nac_width = aircraft.electrofan_nacelle.width
        nac_x_ext = aircraft.electrofan_nacelle.x_ext
        nac_y_ext = aircraft.electrofan_nacelle.y_ext
        nac_z_ext = aircraft.electrofan_nacelle.z_ext
        if (aircraft.electrofan_engine.n_engine==4):
            nac_x_int = aircraft.electrofan_nacelle.x_int
            nac_y_int = aircraft.electrofan_nacelle.y_int
            nac_z_int = aircraft.electrofan_nacelle.z_int
    elif (aircraft.propulsion.architecture=="PTE1"):
        nac_length = aircraft.turbofan_nacelle.length
        nac_height = aircraft.turbofan_nacelle.width
        nac_width = aircraft.turbofan_nacelle.width
        nac_x_ext = aircraft.turbofan_nacelle.x_ext
        nac_y_ext = aircraft.turbofan_nacelle.y_ext
        nac_z_ext = aircraft.turbofan_nacelle.z_ext
    else:
        raise Exception("propulsion.architecture index is out of range")

    e_nac_length = aircraft.rear_electric_nacelle.length
    e_nac_height = aircraft.rear_electric_nacelle.width
    e_nac_width = aircraft.rear_electric_nacelle.width
    e_nac_x_axe = aircraft.rear_electric_nacelle.x_axe
    e_nac_y_axe = aircraft.rear_electric_nacelle.y_axe
    e_nac_z_axe = aircraft.rear_electric_nacelle.z_axe

    r_nose = 0.15       # Fuselage length ratio of nose evolutive part
    r_cone = 0.35       # Fuselage length ratio of tail cone evolutive part

    nose,nose2,cone,cyl = get_shape()

    cyl_yz = np.stack([cyl[0:,0]*fus_width , cyl[0:,1]*fus_height , cyl[0:,2]*fus_height], axis=1)

    fus_front = np.vstack([np.stack([cyl_yz[0:,0] , cyl_yz[0:,1]],axis=1) , np.stack([cyl_yz[::-1,0] , cyl_yz[::-1,2]],axis=1)])

    nose_xz = np.stack([nose[0:,0]*fus_length*r_nose , nose[0:,1]*fus_height , nose[0:,2]*fus_height], axis=1)
    cone_xz = np.stack([(1-r_cone)*fus_length + cone[0:,0]*fus_length*r_cone , cone[0:,1]*fus_height , cone[0:,2]*fus_height], axis=1)
    fus_xz = np.vstack([nose_xz , cone_xz])

    fus_side = np.vstack([np.stack([fus_xz[0:-2,0] , fus_xz[0:-2,1]],axis=1) , np.stack([fus_xz[:0:-1,0] , fus_xz[:0:-1,2]],axis=1)])

    nose_xy = np.stack([nose[0:,0]*fus_length*r_nose , nose[0:,3]*fus_width , nose[0:,4]*fus_width], axis=1)
    cone_xy = np.stack([(1-r_cone)*fus_length + cone[0:,0]*fus_length*r_cone , cone[0:,3]*fus_width , cone[0:,4]*fus_width], axis=1)
    fus_xy = np.vstack([nose_xy , cone_xy])

    fus_top = np.vstack([np.stack([fus_xy[1:-2,0]  , fus_xy[1:-2,1]],axis=1) , np.stack([fus_xy[:0:-1,0] , fus_xy[:0:-1,2]],axis=1)])
    # HTP shape
    #-----------------------------------------------------------------------------------------------------------

    htp_xy = np.array([[htp_x_axe           ,  0            ] ,
                       [htp_x_tip           ,  0.5*htp_span ] ,
                       [htp_x_tip+htp_c_tip ,  0.5*htp_span ] ,
                       [htp_x_axe+htp_c_axe ,  0            ] ,
                       [htp_x_tip+htp_c_tip , -0.5*htp_span ] ,
                       [htp_x_tip           , -0.5*htp_span ] ,
                       [htp_x_axe           ,  0            ]])

    htp_xz = np.array([[htp_x_tip              , htp_z_tip                                ] ,
                       [htp_x_tip+0.1*htp_c_tip , htp_z_tip+0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip+0.7*htp_c_tip , htp_z_tip+0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip+htp_c_tip     , htp_z_tip                               ] ,
                       [htp_x_tip+0.7*htp_c_tip , htp_z_tip-0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip+0.1*htp_c_tip , htp_z_tip-0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip               , htp_z_tip                               ] ,
                       [htp_x_axe               , htp_z_axe                               ] ,
                       [htp_x_axe+0.1*htp_c_axe , htp_z_axe-0.5*htp_t_o_c*htp_c_axe ] ,
                       [htp_x_axe+0.7*htp_c_axe , htp_z_axe-0.5*htp_t_o_c*htp_c_axe ] ,
                       [htp_x_axe+htp_c_axe     , htp_z_axe                               ] ,
                       [htp_x_tip+htp_c_tip     , htp_z_tip                               ] ,
                       [htp_x_tip+0.7*htp_c_tip , htp_z_tip-0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip+0.1*htp_c_tip , htp_z_tip-0.5*htp_t_o_c*htp_c_tip ] ,
                       [htp_x_tip               , htp_z_tip                               ]])

    htp_yz = np.array([[ 0           , htp_z_axe-0.5*fus_height                                                           ] ,
                       [ 0.5*htp_span , htp_z_axe-0.5*fus_height+0.5*htp_span*np.tan(htp_dihedral)                           ] ,
                       [ 0.5*htp_span , htp_z_axe-0.5*fus_height+0.5*htp_span*np.tan(htp_dihedral)-htp_t_o_c*htp_c_tip ] ,
                       [ 0            , htp_z_axe-0.5*fus_height-htp_t_o_c*htp_c_axe                                ] ,
                       [-0.5*htp_span , htp_z_axe-0.5*fus_height+0.5*htp_span*np.tan(htp_dihedral)-htp_t_o_c*htp_c_tip ] ,
                       [-0.5*htp_span , htp_z_axe-0.5*fus_height+0.5*htp_span*np.tan(htp_dihedral)                           ] ,
                       [ 0            , htp_z_axe-0.5*fus_height                                                          ]])

    # VTP shape
    #-----------------------------------------------------------------------------------------------------------

    vtp_xz = np.array([[vtp_x_root           , vtp_z_root  ] ,
                       [vtp_x_tip             , vtp_z_tip  ] ,
                       [vtp_x_tip+vtp_c_tip   , vtp_z_tip  ] ,
                       [vtp_x_root+vtp_c_root , vtp_z_root ]])

    vtp_xy = np.array([[vtp_x_root               ,  0                               ] ,
                       [vtp_x_root+0.1*vtp_c_root ,  0.5*vtp_t_o_c*vtp_c_root ] ,
                       [vtp_x_root+0.7*vtp_c_root ,  0.5*vtp_t_o_c*vtp_c_root ] ,
                       [vtp_x_root+vtp_c_root     ,  0                              ] ,
                       [vtp_x_root+0.7*vtp_c_root , -0.5*vtp_t_o_c*vtp_c_root ] ,
                       [vtp_x_root+0.1*vtp_c_root , -0.5*vtp_t_o_c*vtp_c_root ] ,
                       [vtp_x_root                ,  0                              ] ,
                       [vtp_x_tip                 ,  0                              ] ,
                       [vtp_x_tip+0.1*vtp_c_tip   ,  0.5*vtp_t_o_c*vtp_c_tip  ] ,
                       [vtp_x_tip+0.7*vtp_c_tip   ,  0.5*vtp_t_o_c*vtp_c_tip  ] ,
                       [vtp_x_tip+vtp_c_tip       ,  0                              ] ,
                       [vtp_x_tip+0.7*vtp_c_tip   , -0.5*vtp_t_o_c*vtp_c_tip  ] ,
                       [vtp_x_tip+0.1*vtp_c_tip   , -0.5*vtp_t_o_c*vtp_c_tip  ] ,
                       [vtp_x_tip                 ,  0                              ]])


    vtp_yz = np.array([[ 0.5*vtp_t_o_c*vtp_c_root , -0.5*fus_height+vtp_z_root ],
                       [ 0.5*vtp_t_o_c*vtp_c_tip  , -0.5*fus_height+vtp_z_tip  ],
                       [-0.5*vtp_t_o_c*vtp_c_tip  , -0.5*fus_height+vtp_z_tip  ],
                       [-0.5*vtp_t_o_c*vtp_c_root , -0.5*fus_height+vtp_z_root ]])

    # wing_ shape
    #-----------------------------------------------------------------------------------------------------------

    wing_xy = np.array([[wing_x_root             ,  wing_y_root  ] ,
                        [wing_x_tip              ,  wing_y_tip  ] ,
                        [wing_x_tip+wing_c_tip   ,  wing_y_tip  ] ,
                        [wing_x_kink+wing_c_kink ,  wing_y_kink ] ,
                        [wing_x_root+wing_c_root ,  wing_y_root  ] ,
                        [wing_x_root+wing_c_root , -wing_y_root  ] ,
                        [wing_x_kink+wing_c_kink , -wing_y_kink ] ,
                        [wing_x_tip+wing_c_tip   , -wing_y_tip  ] ,
                        [wing_x_tip              , -wing_y_tip  ] ,
                        [wing_x_root             , -wing_y_root  ] ,
                        [wing_x_root             ,  wing_y_root  ]])

    wing_yz = np.array([[ wing_y_root  , -0.5*fus_height+wing_z_root                        ] ,
                        [ wing_y_kink , -0.5*fus_height+wing_z_kink                        ] ,
                        [ wing_y_tip  , -0.5*fus_height+wing_z_tip                         ] ,
                        [ wing_y_tip  , -0.5*fus_height+wing_z_tip+wing_toc_t*wing_c_tip   ] ,
                        [ wing_y_kink , -0.5*fus_height+wing_z_kink+wing_toc_k*wing_c_kink ] ,
                        [ wing_y_root  , -0.5*fus_height+wing_z_root+wing_toc_r*wing_c_root ] ,
                        [-wing_y_root  , -0.5*fus_height+wing_z_root+wing_toc_r*wing_c_root ] ,
                        [-wing_y_kink , -0.5*fus_height+wing_z_kink+wing_toc_k*wing_c_kink ] ,
                        [-wing_y_tip  , -0.5*fus_height+wing_z_tip+wing_toc_t*wing_c_tip   ] ,
                        [-wing_y_tip  , -0.5*fus_height+wing_z_tip                         ] ,
                        [-wing_y_kink , -0.5*fus_height+wing_z_kink                        ] ,
                        [-wing_y_root  , -0.5*fus_height+wing_z_root                        ] ,
                        [ wing_y_root  , -0.5*fus_height+wing_z_root                        ]])

    wing_xz = np.array([[wing_x_tip                  , wing_z_tip+wing_toc_t*wing_c_tip                           ] ,
                        [wing_x_tip+0.1*wing_c_tip   , wing_z_tip+wing_toc_t*wing_c_tip-0.5*wing_toc_t*wing_c_tip ] ,
                        [wing_x_tip+0.7*wing_c_tip   , wing_z_tip+wing_toc_t*wing_c_tip-0.5*wing_toc_t*wing_c_tip ] ,
                        [wing_x_tip+wing_c_tip       , wing_z_tip+wing_toc_t*wing_c_tip                           ] ,
                        [wing_x_tip+0.7*wing_c_tip   , wing_z_tip+wing_toc_t*wing_c_tip+0.5*wing_toc_t*wing_c_tip ] ,
                        [wing_x_tip+0.1*wing_c_tip   , wing_z_tip+wing_toc_t*wing_c_tip+0.5*wing_toc_t*wing_c_tip ] ,
                        [wing_x_tip                  , wing_z_tip+wing_toc_t*wing_c_tip                           ] ,
                        [wing_x_kink                 , wing_z_kink+0.5*wing_toc_k*wing_c_kink                     ] ,
                        [wing_x_root                 , wing_z_root+0.5*wing_toc_r*wing_c_root                     ] ,
                        [wing_x_root+0.1*wing_c_root , wing_z_root                                                ] ,
                        [wing_x_root+0.7*wing_c_root , wing_z_root                                                ] ,
                        [wing_x_root+wing_c_root     , wing_z_root+0.5*wing_toc_r*wing_c_root                     ] ,
                        [wing_x_kink+wing_c_kink     , wing_z_kink+0.5*wing_toc_k*wing_c_kink                     ] ,
                        [wing_x_tip+wing_c_tip       , wing_z_tip+wing_toc_t*wing_c_tip                           ]])


    # External engine shape
    #-----------------------------------------------------------------------------------------------------------
    nac_xz_ext,nac_xy_ext,nac_yz_ext,fan_yz_ext = nacelle_shape(0.5*fus_height, \
                           nac_x_ext,nac_y_ext,nac_z_ext,nac_width,nac_height,nac_length,cyl)

    if (aircraft.propulsion.architecture=="TF"):
        if (aircraft.turbofan_engine.n_engine==4):
            nac_xz_int,nac_xy_int,nac_yz_int,fan_yz_int = nacelle_shape(0.5*fus_height, \
                           nac_x_int,nac_y_int,nac_z_int,nac_width,nac_height,nac_length,cyl)

    # e-Pods
    #-----------------------------------------------------------------------------------------------------------
    if (aircraft.propulsion.architecture=="PTE1"):
        e_nac_xz = np.array([[e_nac_x_axe                 , e_nac_z_axe+0.5*fus_height+0.4*e_nac_width ] ,
                            [e_nac_x_axe+0.1*e_nac_length , e_nac_z_axe+0.5*fus_height+0.5*e_nac_width ] ,
                            [e_nac_x_axe+0.7*e_nac_length , e_nac_z_axe+0.5*fus_height+0.5*e_nac_width ] ,
                            [e_nac_x_axe+e_nac_length     , e_nac_z_axe+0.5*fus_height+0.4*e_nac_width ] ,
                            [e_nac_x_axe+e_nac_length     , e_nac_z_axe+0.5*fus_height-0.4*e_nac_width ] ,
                            [e_nac_x_axe+0.7*e_nac_length , e_nac_z_axe+0.5*fus_height-0.5*e_nac_width ] ,
                            [e_nac_x_axe+0.1*e_nac_length , e_nac_z_axe+0.5*fus_height-0.5*e_nac_width ] ,
                            [e_nac_x_axe                 , e_nac_z_axe+0.5*fus_height-0.4*e_nac_width ] ,
                            [e_nac_x_axe                 , e_nac_z_axe+0.5*fus_height+0.4*e_nac_width ]])

        e_nac_xy = np.array([[e_nac_x_axe                 ,  0.4*e_nac_width ] ,
                            [e_nac_x_axe+0.1*e_nac_length ,  0.5*e_nac_width ] ,
                            [e_nac_x_axe+0.7*e_nac_length ,  0.5*e_nac_width ] ,
                            [e_nac_x_axe+e_nac_length     ,  0.4*e_nac_width ] ,
                            [e_nac_x_axe+e_nac_length     , -0.4*e_nac_width ] ,
                            [e_nac_x_axe+0.7*e_nac_length , -0.5*e_nac_width ] ,
                            [e_nac_x_axe+0.1*e_nac_length , -0.5*e_nac_width ] ,
                            [e_nac_x_axe                 , -0.4*e_nac_width ] ,
                            [e_nac_x_axe                 ,  0.4*e_nac_width ]])

        e_d_nac_yz = np.stack([cyl[0:,0]*e_nac_width , cyl[0:,1]*e_nac_width , cyl[0:,2]*e_nac_width], axis=1)

        e_d_fan_yz = np.stack([cyl[0:,0]*0.80*e_nac_width , cyl[0:,1]*0.80*e_nac_width , cyl[0:,2]*0.80*e_nac_width], axis=1)

        e_nac_yz = np.vstack([np.stack([e_nac_y_axe+e_d_nac_yz[0:,0] , e_nac_z_axe+e_d_nac_yz[0:,1]],axis=1) ,
                                 np.stack([e_nac_y_axe+e_d_nac_yz[::-1,0] , e_nac_z_axe+e_d_nac_yz[::-1,2]],axis=1)])

        e_fan_yz = np.vstack([np.stack([e_nac_y_axe+e_d_fan_yz[0:,0] , e_nac_z_axe+e_d_fan_yz[0:,1]],axis=1) ,
                                 np.stack([e_nac_y_axe+e_d_fan_yz[::-1,0] , e_nac_z_axe+e_d_fan_yz[::-1,2]],axis=1)])

    # Drawing_ box
    #-----------------------------------------------------------------------------------------------------------
    fig,axes = plt.subplots(1,1)
    fig.canvas.set_window_title(window_title)
    fig.suptitle(plot_title, fontsize=14)
    axes.set_aspect('equal', 'box')
    plt.plot(np.array([0,100,100,0,0]), np.array([0,0,100,100,0]))      # Draw a square box of 100m side

    xTopView = 50 - (aircraft.wing.x_mac + 0.25*aircraft.wing.mac)      # Top view positionning
    yTopView = 50

    xSideView = 50 - (aircraft.wing.x_mac+0.25*aircraft.wing.mac)       # Top view positionning
    ySideView = 82

    xFrontView = 50
    yFrontView = 10

    # Draw top view
    #-----------------------------------------------------------------------------------------------------------
    if (aircraft.propulsion.architecture!="TP"):
        plt.plot(xTopView+nac_xy_ext[0:,0], yTopView+nac_xy_ext[0:,1], color="grey", zorder=3)        # Left nacelle top view
        plt.plot(xTopView+nac_xy_ext[0:,0], yTopView-nac_xy_ext[0:,1], color="grey", zorder=3)        # Right nacelle top view
        if (aircraft.turbofan_engine.n_engine==4):
            plt.plot(xTopView+nac_xy_int[0:,0], yTopView+nac_xy_int[0:,1], color="grey", zorder=3)        # Left nacelle top view
            plt.plot(xTopView+nac_xy_int[0:,0], yTopView-nac_xy_int[0:,1], color="grey", zorder=3)        # Right nacelle top view
        plt.fill(xTopView+wing_xy[0:,0], yTopView+wing_xy[0:,1], color="white", zorder=4)     # wing_ top view
        plt.plot(xTopView+wing_xy[0:,0], yTopView+wing_xy[0:,1], color="grey", zorder=4)      # wing_ top view
    elif (aircraft.propulsion.architecture=="TP"):
        plt.fill(xTopView+wing_xy[0:,0], yTopView+wing_xy[0:,1], color="white", zorder=3)     # wing_ top view
        plt.plot(xTopView+wing_xy[0:,0], yTopView+wing_xy[0:,1], color="grey", zorder=3)      # wing_ top view
        plt.fill(xTopView+nac_xy_ext[0:,0], yTopView+nac_xy_ext[0:,1], color="white", zorder=4)        # Left nacelle top view
        plt.plot(xTopView+nac_xy_ext[0:,0], yTopView+nac_xy_ext[0:,1], color="grey", zorder=4)        # Left nacelle top view
        plt.fill(xTopView+nac_xy_ext[0:,0], yTopView-nac_xy_ext[0:,1], color="white", zorder=4)        # Right nacelle top view
        plt.plot(xTopView+nac_xy_ext[0:,0], yTopView-nac_xy_ext[0:,1], color="grey", zorder=4)        # Right nacelle top view

    if (aircraft.horizontal_tail.attachment==1):
        plt.plot(xTopView+htp_xy[0:,0], yTopView+htp_xy[0:,1], "grey", zorder=1)      # htp_ top view (Classic or Vtail)

    if (aircraft.wing.attachment==1):
        plt.fill(xTopView+fus_top[0:,0], yTopView+fus_top[0:,1], color="white", zorder=6)   # fuselage top view
        plt.plot(xTopView+fus_top[0:,0], yTopView+fus_top[0:,1], "grey", zorder=6)          # fuselage top view
    elif (aircraft.wing.attachment==2):
        plt.fill(xTopView+fus_top[0:,0], yTopView+fus_top[0:,1], color="white", zorder=2)   # fuselage top view
        plt.plot(xTopView+fus_top[0:,0], yTopView+fus_top[0:,1], "grey", zorder=2)          # fuselage top view

    plt.plot(xTopView+vtp_xy[0:,0], yTopView+vtp_xy[0:,1], "grey", zorder=8)            # vtp top view

    if (aircraft.horizontal_tail.attachment==2):
        plt.plot(xTopView+htp_xy[0:,0], yTopView+htp_xy[0:,1], "grey", zorder=9)      # htp_ top view (T-tail)

    if (aircraft.propulsion.architecture=="PTE1"):
        plt.plot(xTopView+e_nac_xy[0:,0], yTopView+e_nac_xy[0:,1], color="grey", zorder=7)        # r_nacelle top view

    # Draw side view
    #-----------------------------------------------------------------------------------------------------------
    plt.plot(xSideView+vtp_xz[0:,0], ySideView+vtp_xz[0:,1], color="grey", zorder=1)      # vtp_ side view

    plt.fill(xSideView+fus_side[0:,0], ySideView+fus_side[0:,1], color="white", zorder=2) # fuselage side view
    plt.plot(xSideView+fus_side[0:,0], ySideView+fus_side[0:,1], color="grey", zorder=3)  # fuselage side view

    if (aircraft.propulsion.architecture=="PTE1"):
        plt.fill(xSideView+e_nac_xz[0:,0], ySideView+e_nac_xz[0:,1], color="white", zorder=4)   # r_nacelle side view
        plt.plot(xSideView+e_nac_xz[0:,0], ySideView+e_nac_xz[0:,1], color="grey", zorder=5)    # r_nacelle side view

    if (aircraft.propulsion.architecture!="TP"):
        plt.fill(xSideView+wing_xz[0:,0], ySideView+wing_xz[0:,1], color="white", zorder=4)   # wing_ side view
        plt.plot(xSideView+wing_xz[0:,0], ySideView+wing_xz[0:,1], color="grey", zorder=5)    # wing_ side view
        if (aircraft.turbofan_engine.n_engine==4):
            plt.fill(xSideView+nac_xz_int[0:,0], ySideView+nac_xz_int[0:,1], color="white", zorder=6)     # nacelle side view
            plt.plot(xSideView+nac_xz_int[0:,0], ySideView+nac_xz_int[0:,1], color="grey", zorder=6)      # nacelle side view
        plt.fill(xSideView+nac_xz_ext[0:,0], ySideView+nac_xz_ext[0:,1], color="white", zorder=7)     # nacelle side view
        plt.plot(xSideView+nac_xz_ext[0:,0], ySideView+nac_xz_ext[0:,1], color="grey", zorder=7)      # nacelle side view
    elif (aircraft.propulsion.architecture=="TP"):
        plt.fill(xSideView+nac_xz_ext[0:,0], ySideView+nac_xz_ext[0:,1], color="white", zorder=4)     # nacelle side view
        plt.plot(xSideView+nac_xz_ext[0:,0], ySideView+nac_xz_ext[0:,1], color="grey", zorder=5)      # nacelle side view
        plt.fill(xSideView+wing_xz[0:,0], ySideView+wing_xz[0:,1], color="white", zorder=6)   # wing_ side view
        plt.plot(xSideView+wing_xz[0:,0], ySideView+wing_xz[0:,1], color="grey", zorder=7)    # wing_ side view

    plt.fill(xSideView+htp_xz[0:,0], ySideView+htp_xz[0:,1], color="white", zorder=4)     # htp_ side view
    plt.plot(xSideView+htp_xz[0:,0], ySideView+htp_xz[0:,1], color="grey", zorder=5)      # htp_ side view

    # Draw front view
    #-----------------------------------------------------------------------------------------------------------
    plt.plot(xFrontView+vtp_yz[0:,0], yFrontView+vtp_yz[0:,1], color="grey", zorder=1)     # vtp_ front view
    plt.plot(xFrontView+htp_yz[0:,0], yFrontView+htp_yz[0:,1], color="grey", zorder=1)     # htp_ front view

    plt.plot(xFrontView+wing_yz[0:,0], yFrontView+wing_yz[0:,1], color="grey", zorder=2)   # wing_ front view

    if (aircraft.propulsion.architecture=="PTE1"):
        plt.plot(xFrontView+e_nac_yz[0:,0], yFrontView+e_nac_yz[0:,1], color="grey", zorder=3)    # r_nacelle front view
        plt.plot(xFrontView+e_fan_yz[0:,0], yFrontView+e_fan_yz[0:,1], color="grey", zorder=3)    # e inlet front view

    plt.fill(xFrontView+fus_front[0:,0], yFrontView+fus_front[0:,1], color="white", zorder=4)   # fuselage front view
    plt.plot(xFrontView+fus_front[0:,0], yFrontView+fus_front[0:,1], color="grey", zorder=5)    # fuselage front view

    plt.fill(xFrontView+nac_yz_ext[0:,0], yFrontView+nac_yz_ext[0:,1], color="white", zorder=6)   # Left nacelle front view
    plt.plot(xFrontView+nac_yz_ext[0:,0], yFrontView+nac_yz_ext[0:,1], color="grey", zorder=7)    # Left nacelle front view
    plt.plot(xFrontView+fan_yz_ext[0:,0], yFrontView+fan_yz_ext[0:,1], color="grey", zorder=8)    # Left Inlet front view
    if (aircraft.turbofan_engine.n_engine==0):
        plt.fill(xFrontView+nac_yz_int[0:,0], yFrontView+nac_yz_int[0:,1], color="white", zorder=6)   # Left nacelle front view
        plt.plot(xFrontView+nac_yz_int[0:,0], yFrontView+nac_yz_int[0:,1], color="grey", zorder=7)    # Left nacelle front view
        plt.plot(xFrontView+fan_yz_int[0:,0], yFrontView+fan_yz_int[0:,1], color="grey", zorder=8)    # Left Inlet front view

    plt.fill(xFrontView-nac_yz_ext[0:,0], yFrontView+nac_yz_ext[0:,1], color="white", zorder=6)   # Right nacelle front view
    plt.plot(xFrontView-nac_yz_ext[0:,0], yFrontView+nac_yz_ext[0:,1], color="grey", zorder=7)    # Right nacelle front view
    plt.plot(xFrontView-fan_yz_ext[0:,0], yFrontView+fan_yz_ext[0:,1], color="grey", zorder=8)    # Right Inlet front view
    if (aircraft.turbofan_engine.n_engine==4):
        plt.fill(xFrontView-nac_yz_int[0:,0], yFrontView+nac_yz_int[0:,1], color="white", zorder=6)   # Right nacelle front view
        plt.plot(xFrontView-nac_yz_int[0:,0], yFrontView+nac_yz_int[0:,1], color="grey", zorder=7)    # Right nacelle front view
        plt.plot(xFrontView-fan_yz_int[0:,0], yFrontView+fan_yz_int[0:,1], color="grey", zorder=8)    # Right Inlet front view

    plt.show()

    return


#===========================================================================================================
def nacelle_shape(z_shift,nac_x,nac_y,nac_z,nac_width,nac_height,nac_length,cyl):

    nac_xz = np.array([[nac_x                , z_shift+nac_z+0.4*nac_height ] ,
                       [nac_x+0.1*nac_length , z_shift+nac_z+0.5*nac_height ] ,
                       [nac_x+0.5*nac_length , z_shift+nac_z+0.5*nac_height ] ,
                       [nac_x+nac_length     , z_shift+nac_z+0.3*nac_height ] ,
                       [nac_x+nac_length     , z_shift+nac_z-0.3*nac_height ] ,
                       [nac_x+0.5*nac_length , z_shift+nac_z-0.5*nac_height ] ,
                       [nac_x+0.1*nac_length , z_shift+nac_z-0.5*nac_height ] ,
                       [nac_x                , z_shift+nac_z-0.4*nac_height ] ,
                       [nac_x                , z_shift+nac_z+0.4*nac_height ]])
    
    nac_xy = np.array([[nac_x                , nac_y+0.4*nac_width ] ,
                       [nac_x+0.1*nac_length , nac_y+0.5*nac_width ] ,
                       [nac_x+0.5*nac_length , nac_y+0.5*nac_width ] ,
                       [nac_x+nac_length     , nac_y+0.3*nac_width ] ,
                       [nac_x+nac_length     , nac_y-0.3*nac_width ] ,
                       [nac_x+0.5*nac_length , nac_y-0.5*nac_width ] ,
                       [nac_x+0.1*nac_length , nac_y-0.5*nac_width ] ,
                       [nac_x                , nac_y-0.4*nac_width ] ,
                       [nac_x                , nac_y+0.4*nac_width ]])
    
    d_nac_yz = np.stack([cyl[0:,0]*nac_width , cyl[0:,1]*nac_height , cyl[0:,2]*nac_height], axis=1)
    
    d_fan_yz = np.stack([cyl[0:,0]*0.80*nac_width , cyl[0:,1]*0.80*nac_height , cyl[0:,2]*0.80*nac_height], axis=1)
    
    nac_yz = np.vstack([np.stack([nac_y+d_nac_yz[0:,0] , nac_z+d_nac_yz[0:,1]],axis=1) ,
                           np.stack([nac_y+d_nac_yz[::-1,0] , nac_z+d_nac_yz[::-1,2]],axis=1)])
    
    fan_yz = np.vstack([np.stack([nac_y+d_fan_yz[0:,0] , nac_z+d_fan_yz[0:,1]],axis=1) ,
                           np.stack([nac_y+d_fan_yz[::-1,0] , nac_z+d_fan_yz[::-1,2]],axis=1)])

    return nac_xz,nac_xy,nac_yz,fan_yz

#===========================================================================================================
def get_shape():
    """
    Total aircraft drag with the assumption that the wing takes all the lift
    """

    nose1 = np.array([[ 0.0000 , 0.4453 , 0.4453 , 0.0000 ,  0.0000 ] ,
                      [ 0.0050 , 0.4733 , 0.4112 , 0.0335 , -0.0335 ] ,
                      [ 0.0191 , 0.5098 , 0.3833 , 0.0646 , -0.0646 ] ,
                      [ 0.0624 , 0.5718 , 0.3188 , 0.1196 , -0.1196 ] ,
                      [ 0.1355 , 0.6278 , 0.2531 , 0.1878 , -0.1878 ] ,
                      [ 0.1922 , 0.7263 , 0.2142 , 0.2297 , -0.2297 ] ,
                      [ 0.2773 , 0.8127 , 0.1631 , 0.2859 , -0.2859 ] ,
                      [ 0.4191 , 0.8906 , 0.0962 , 0.3624 , -0.3624 ] ,
                      [ 0.5610 , 0.9392 , 0.0536 , 0.4211 , -0.4211 ] ,
                      [ 0.7738 , 0.9818 , 0.0122 , 0.4761 , -0.4761 ] ,
                      [ 0.9156 , 0.9976 , 0.0025 , 0.4976 , -0.4976 ] ,
                      [ 1.0000 , 1.0000 , 0.0000 , 0.5000 , -0.5000 ]])

    nose2 = np.array([[ 0.0000 , 0.3339 , 0.3339 , 0.0000 ,  0.0000 ] ,
                      [ 0.0050 , 0.3848 , 0.3084 , 0.0335 , -0.0335 ] ,
                      [ 0.0150 , 0.4253 , 0.2881 , 0.0652 , -0.0652 ] ,
                      [ 0.0500 , 0.5033 , 0.2490 , 0.1101 , -0.1101 ] ,
                      [ 0.1000 , 0.5811 , 0.2100 , 0.1585 , -0.1585 ] ,
                      [ 0.1800 , 0.6808 , 0.1600 , 0.2215 , -0.2215 ] ,
                      [ 0.2773 , 0.7704 , 0.1151 , 0.2859 , -0.2859 ] ,
                      [ 0.4191 , 0.8562 , 0.0721 , 0.3624 , -0.3624 ] ,
                      [ 0.5610 , 0.9198 , 0.0402 , 0.4211 , -0.4211 ] ,
                      [ 0.7738 , 0.9816 , 0.0092 , 0.4761 , -0.4761 ] ,
                      [ 0.9156 , 0.9962 , 0.0019 , 0.4976 , -0.4976 ] ,
                      [ 1.0000 , 1.0000 , 0.0000 , 0.5000 , -0.5000 ]])

    cone = np.array([[ 0.0000 , 1.0000 , 0.0000 , 0.5000 , -0.5000 ] ,
                     [ 0.0213 , 1.0000 , 0.0082 , 0.5029 , -0.5029 ] ,
                     [ 0.0638 , 1.0000 , 0.0230 , 0.4956 , -0.4956 ] ,
                     [ 0.1064 , 1.0000 , 0.0393 , 0.4875 , -0.4875 ] ,
                     [ 0.1489 , 1.0000 , 0.0556 , 0.4794 , -0.4794 ] ,
                     [ 0.1915 , 1.0000 , 0.0786 , 0.4720 , -0.4720 ] ,
                     [ 0.2766 , 1.0000 , 0.1334 , 0.4566 , -0.4566 ] ,
                     [ 0.3617 , 1.0000 , 0.1964 , 0.4330 , -0.4330 ] ,
                     [ 0.4894 , 1.0000 , 0.3024 , 0.3822 , -0.3822 ] ,
                     [ 0.6170 , 1.0000 , 0.4159 , 0.3240 , -0.3240 ] ,
                     [ 0.7447 , 1.0000 , 0.5374 , 0.2577 , -0.2577 ] ,
                     [ 0.8723 , 1.0000 , 0.6627 , 0.1834 , -0.1834 ] ,
                     [ 0.8936 , 0.9963 , 0.6901 , 0.1679 , -0.1679 ] ,
                     [ 0.9149 , 0.9881 , 0.7139 , 0.1524 , -0.1524 ] ,
                     [ 0.9362 , 0.9800 , 0.7413 , 0.1333 , -0.1333 ] ,
                     [ 0.9574 , 0.9652 , 0.7687 , 0.1097 , -0.1097 ] ,
                     [ 0.9787 , 0.9533 , 0.8043 , 0.0788 , -0.0788 ] ,
                     [ 0.9894 , 0.9377 , 0.8280 , 0.0589 , -0.0589 ] ,
                     [ 1.0000 , 0.9103 , 0.8784 , 0.0162 , -0.0162 ]])

    cyl = np.array([[  0.5000000 , 0.0000000 ,  0.0000000 ] ,
                    [  0.4903926 , 0.0975452 , -0.0975452 ] ,
                    [  0.4619398 , 0.1913417 , -0.1913417 ] ,
                    [  0.4157348 , 0.2777851 , -0.2777851 ] ,
                    [  0.3535534 , 0.3535534 , -0.3535534 ] ,
                    [  0.2777851 , 0.4157348 , -0.4157348 ] ,
                    [  0.1913417 , 0.4619398 , -0.4619398 ] ,
                    [  0.0975452 , 0.4903926 , -0.4903926 ] ,
                    [  0.0000000 , 0.5000000 , -0.5000000 ] ,
                    [- 0.0975452 , 0.4903926 , -0.4903926 ] ,
                    [- 0.1913417 , 0.4619398 , -0.4619398 ] ,
                    [- 0.2777851 , 0.4157348 , -0.4157348 ] ,
                    [- 0.3535534 , 0.3535534 , -0.3535534 ] ,
                    [- 0.4157348 , 0.2777851 , -0.2777851 ] ,
                    [- 0.4619398 , 0.1913417 , -0.1913417 ] ,
                    [- 0.4903926 , 0.0975452 , -0.0975452 ] ,
                    [- 0.5000000 , 0.0000000 ,  0.0000000 ]])

    return nose1,nose2,cone,cyl

