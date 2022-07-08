from conflux.bsg.Constants import *

def screening_potential(Z):
    """Get shift in electric potential at nuclear radius due to atomic electrons

    :param Z: Proton number of the nuclear state
    """
    
    if Z < 1:
        return 0

    a, b = screening_parameters(Z)

    return 2*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

def screening_parameters(Z):
    """ Returns the fit parameters of the atomic potential
    by Salvat et al., Physical Review A 36 (1987) 467-474
    for Z <= 92. If Z is larger, the Moliere potential is used

    :param Z: the proton number of the atom
    """
    a = []
    b = []
    if Z < 0:
        Z = -Z  # for Z<0 if beta + transition

    if Z > 92:  # Moliere's potential
        c = 0.88534 * pow(Z * 1., -1. / 3.)
        a.append(0.1)
        a.append(0.55)
        a.append(0.35)
        b.append(6.0 / c)
        b.append(1.2 / c)
        b.append(0.3 / c)

    else:
        # Screening for H1
        if Z == 1:
            a.append(-184.39)
            a.append(185.39)
            a.append(1. - (a[0] + a[1]))
            b.append(2.0027)
            b.append(1.9973)
            b.append(0.0000)

        # Screening for He2
        if Z == 2:
            a.append(-0.2259)
            a.append(1.2259)
            a.append(1. - (a[0] + a[1]))
            b.append(5.5272)
            b.append(2.3992)
            b.append(0.0000)

        # Screening for Li3
        if Z == 3:
            a.append(0.6045)
            a.append(0.3955)
            a.append(1. - (a[0] + a[1]))
            b.append(2.8174)
            b.append(0.6625)
            b.append(0.0000)

        # Screening for Be4
        if Z == 4:
            a.append(0.3278)
            a.append(0.6722)
            a.append(1. - (a[0] + a[1]))
            b.append(04.5430)
            b.append(0.9852)
            b.append(0.0000)

        # Screening for B5
        if Z == 5:
            a.append(0.2327)
            a.append(0.7673)
            a.append(1. - (a[0] + a[1]))
            b.append(5.9900)
            b.append(1.2135)
            b.append(0.0000)

        # Screening for C6
        if Z == 6:
            a.append(0.1537)
            a.append(0.8463)
            a.append(1. - (a[0] + a[1]))
            b.append(8.0404)
            b.append(1.4913)
            b.append(0.0000)

        # Screening for N7
        if Z == 7:
            a.append(0.0996)
            a.append(0.9004)
            a.append(1. - (a[0] + a[1]))
            b.append(10.812)
            b.append(1.7687)
            b.append(0.0000)

        # Screening for O8
        if Z == 8:
            a.append(0.0625)
            a.append(0.9375)
            a.append(1. - (a[0] + a[1]))
            b.append(14.823)
            b.append(2.0403)
            b.append(0.0000)

        # Screening for F9
        if Z == 9:
            a.append(0.0368)
            a.append(0.9632)
            a.append(1. - (a[0] + a[1]))
            b.append(21.400)
            b.append(2.3060)
            b.append(0.0000)

        # Screening for Ne10
        if Z == 10:
            a.append(0.0188)
            a.append(0.9812)
            a.append(1 - (a[0] + a[1]))
            b.append(34.999)
            b.append(2.5662)
            b.append(0.0000)

        # Screening for Na11
        if Z == 11:
            a.append(0.7444)
            a.append(0.2556)
            a.append(1. - (a[0] + a[1]))
            b.append(4.1205)
            b.append(0.8718)
            b.append(0.0000)

        # Screening for Mg12
        if Z == 12:
            a.append(0.6423)
            a.append(0.3577)
            a.append(1. - (a[0] + a[1]))
            b.append(4.7266)
            b.append(1.0025)
            b.append(0.0000)

        # Screening for Al13
        if Z == 13:
            a.append(0.6002)
            a.append(0.3998)
            a.append(1. - (a[0] + a[1]))
            b.append(5.1405)
            b.append(1.0153)
            b.append(0.0000)

        # Screening for Si14
        if Z == 14:
            a.append(0.5160)
            a.append(0.4840)
            a.append(1. - (a[0] + a[1]))
            b.append(5.8492)
            b.append(1.1732)
            b.append(0.0000)

        # Screening for P15
        if Z == 15:
            a.append(0.4387)
            a.append(0.5613)
            a.append(1. - (a[0] + a[1]))
            b.append(6.6707)
            b.append(1.3410)
            b.append(0.0000)

        # Screening for S16
        if Z == 16:
            a.append(0.5459)
            a.append(-0.5333)
            a.append(1. - (a[0] + a[1]))
            b.append(6.3703)
            b.append(2.5517)
            b.append(1.6753)

        # Screening for Cl17
        if Z == 17:
            a.append(0.7249)
            a.append(-0.7548)
            a.append(1. - (a[0] + a[1]))
            b.append(6.2118)
            b.append(3.3883)
            b.append(1.8596)

        # Screening for Ar18
        if Z == 18:
            a.append(2.1912)
            a.append(-2.2852)
            a.append(1. - (a[0] + a[1]))
            b.append(5.5470)
            b.append(4.5687)
            b.append(2.0446)

        # Screening for K19
        if Z == 19:
            a.append(0.0486)
            a.append(0.7759)
            a.append(1. - (a[0] + a[1]))
            b.append(30.260)
            b.append(3.1243)
            b.append(0.7326)

        # Screening for Ca20
        if Z == 20:
            a.append(0.5800)
            a.append(0.4200)
            a.append(1. - (a[0] + a[1]))
            b.append(6.3218)
            b.append(1.0094)
            b.append(0.0000)

        # Screening for Sc21
        if Z == 21:
            a.append(0.5543)
            a.append(0.4457)
            a.append(1. - (a[0] + a[1]))
            b.append(6.6328)
            b.append(1.1023)
            b.append(0.0000)

        # Screening for Ti22
        if Z == 22:
            a.append(0.0112)
            a.append(0.6832)
            a.append(1. - (a[0] + a[1]))
            b.append(99.757)
            b.append(4.1286)
            b.append(1.0090)

        # Screening for V23
        if Z == 23:
            a.append(0.0318)
            a.append(0.6753)
            a.append(1. - (a[0] + a[1]))
            b.append(42.533)
            b.append(3.9404)
            b.append(1.0533)

        # Screening for Cr24
        if Z == 24:
            a.append(0.1075)
            a.append(0.7162)
            a.append(1. - (a[0] + a[1]))
            b.append(18.959)
            b.append(3.0638)
            b.append(1.0014)

        # Screening for Mn25
        if Z == 25:
            a.append(0.0498)
            a.append(0.6866)
            a.append(1. - (a[0] + a[1]))
            b.append(31.864)
            b.append(3.7811)
            b.append(1.1279)

        # Screening for Fe26
        if Z == 26:
            a.append(0.0512)
            a.append(0.6995)
            a.append(1. - (a[0] + a[1]))
            b.append(31.825)
            b.append(3.7716)
            b.append(1.1606)

        # Screening for Co27
        if Z == 27:
            a.append(0.0500)
            a.append(0.7142)
            a.append(1. - (a[0] + a[1]))
            b.append(32.915)
            b.append(3.7908)
            b.append(1.1915)

        # Screening for Ni28
        if Z == 28:
            a.append(0.0474)
            a.append(0.7294)
            a.append(1. - (a[0] + a[1]))
            b.append(34.758)
            b.append(3.8299)
            b.append(1.2209)

        # Screening for Cu29
        if Z == 29:
            a.append(0.0771)
            a.append(0.7951)
            a.append(1. - (a[0] + a[1]))
            b.append(25.326)
            b.append(3.3928)
            b.append(1.1426)

        # Screening for Zn30
        if Z == 30:
            a.append(0.0400)
            a.append(0.7590)
            a.append(1. - (a[0] + a[1]))
            b.append(40.343)
            b.append(3.9465)
            b.append(1.2759)

        # Screening for Ga31
        if Z == 31:
            a.append(0.1083)
            a.append(0.7489)
            a.append(1. - (a[0] + a[1]))
            b.append(20.192)
            b.append(3.4733)
            b.append(1.0064)

        # Screening for Ge32
        if Z == 32:
            a.append(0.0610)
            a.append(0.7157)
            a.append(1. - (a[0] + a[1]))
            b.append(29.200)
            b.append(4.1252)
            b.append(1.1845)

        # Screening for As33
        if Z == 33:
            a.append(0.0212)
            a.append(0.6709)
            a.append(1. - (a[0] + a[1]))
            b.append(62.487)
            b.append(4.9502)
            b.append(1.3582)

        # Screening for Se34
        if Z == 34:
            a.append(0.4836)
            a.append(0.5164)
            a.append(1. - (a[0] + a[1]))
            b.append(8.7824)
            b.append(1.6967)
            b.append(0.0000)

        # Screening for Br35
        if Z == 35:
            a.append(0.4504)
            a.append(0.5496)
            a.append(1. - (a[0] + a[1]))
            b.append(9.3348)
            b.append(1.7900)
            b.append(0.0000)

        # Screening for Kr36
        if Z == 36:
            a.append(0.4190)
            a.append(0.5810)
            a.append(1. - (a[0] + a[1]))
            b.append(09.9142)
            b.append(1.8835)
            b.append(0.0000)

        # Screening for Rb37
        if Z == 37:
            a.append(0.1734)
            a.append(0.7253)
            a.append(1. - (a[0] + a[1]))
            b.append(17.166)
            b.append(3.1103)
            b.append(0.7177)

        # Screening for Sr38
        if Z == 38:
            a.append(0.0336)
            a.append(0.7816)
            a.append(1. - (a[0] + a[1]))
            b.append(55.208)
            b.append(4.2842)
            b.append(0.8578)

        # Screening for Y39
        if Z == 39:
            a.append(0.0689)
            a.append(0.7202)
            a.append(1. - (a[0] + a[1]))
            b.append(31.366)
            b.append(4.2412)
            b.append(0.9472)

        # Screening for Zr40
        if Z == 40:
            a.append(0.1176)
            a.append(0.6581)
            a.append(1. - (a[0] + a[1]))
            b.append(22.054)
            b.append(4.0325)
            b.append(1.0181)

        # Screening for Nb41
        if Z == 41:
            a.append(0.2257)
            a.append(0.5821)
            a.append(1. - (a[0] + a[1]))
            b.append(14.240)
            b.append(2.9702)
            b.append(1.0170)

        # Screening for Mo42
        if Z == 42:
            a.append(0.2693)
            a.append(0.5763)
            a.append(1. - (a[0] + a[1]))
            b.append(14.044)
            b.append(2.8611)
            b.append(1.0591)

        # Screening for Tc43
        if Z == 43:
            a.append(0.2201)
            a.append(0.5618)
            a.append(1. - (a[0] + a[1]))
            b.append(15.918)
            b.append(3.3672)
            b.append(1.1548)

        # Screening for Ru44
        if Z == 44:
            a.append(0.2751)
            a.append(0.5943)
            a.append(1. - (a[0] + a[1]))
            b.append(14.314)
            b.append(2.7370)
            b.append(1.1092)

        # Screening for Rh45
        if Z == 45:
            a.append(0.2711)
            a.append(0.6119)
            a.append(1. - (a[0] + a[1]))
            b.append(14.654)
            b.append(2.7183)
            b.append(1.1234)

        # Screening for Pd46
        if Z == 46:
            a.append(0.2784)
            a.append(0.6067)
            a.append(1. - (a[0] + a[1]))
            b.append(14.645)
            b.append(2.6155)
            b.append(1.4318)

        # Screening for Ag47
        if Z == 47:
            a.append(0.2562)
            a.append(0.6505)
            a.append(1. - (a[0] + a[1]))
            b.append(15.5880)
            b.append(2.7412)
            b.append(1.1408)

        # Screening for Cd48

        if Z == 48:
            a.append(0.2271)
            a.append(0.6155)
            a.append(1. - (a[0] + a[1]))
            b.append(16.914)
            b.append(3.0841)
            b.append(1.2619)

        # Screening for In49
        if Z == 49:
            a.append(0.2492)
            a.append(0.6440)
            a.append(1. - (a[0] + a[1]))
            b.append(16.155)
            b.append(2.8819)
            b.append(0.9942)

        # Screening for Sn50
        if Z == 50:
            a.append(0.2153)
            a.append(0.6115)
            a.append(1. - (a[0] + a[1]))
            b.append(17.7930)
            b.append(3.2937)
            b.append(1.1478)

        # Screening for Sb51
        if Z == 51:
            a.append(0.1806)
            a.append(0.5767)
            a.append(1. - (a[0] + a[1]))
            b.append(19.875)
            b.append(3.8092)
            b.append(1.2829)

        # Screening for Te52
        if Z == 52:
            a.append(0.1308)
            a.append(0.5504)
            a.append(1. - (a[0] + a[1]))
            b.append(24.154)
            b.append(4.6119)
            b.append(1.4195)

        # Screening for I53
        if Z == 53:
            a.append(0.0588)
            a.append(0.5482)
            a.append(1. - (a[0] + a[1]))
            b.append(39.996)
            b.append(5.9132)
            b.append(1.5471)

        # Screening for Xe54
        if Z == 54:
            a.append(0.4451)
            a.append(0.5549)
            a.append(1. - (a[0] + a[1]))
            b.append(11.8050)
            b.append(1.7967)
            b.append(0.0000)

        # Screening for Cs55
        if Z == 55:
            a.append(0.2708)
            a.append(0.6524)
            a.append(1. - (a[0] + a[1]))
            b.append(16.591)
            b.append(2.6964)
            b.append(0.6814)

        # Screening for Ba56
        if Z == 56:
            a.append(0.1728)
            a.append(0.6845)
            a.append(1. - (a[0] + a[1]))
            b.append(22.397)
            b.append(3.4595)
            b.append(0.8073)

        # Screening for La57
        if Z == 57:
            a.append(0.1947)
            a.append(0.6384)
            a.append(1. - (a[0] + a[1]))
            b.append(20.764)
            b.append(3.4657)
            b.append(0.8911)

        # Screening for Ce58
        if Z == 58:
            a.append(0.1913)
            a.append(0.6467)
            a.append(1. - (a[0] + a[1]))
            b.append(21.235)
            b.append(3.4819)
            b.append(0.9011)

        # Screening for Pr59
        if Z == 59:
            a.append(0.1868)
            a.append(0.6558)
            a.append(1. - (a[0] + a[1]))
            b.append(21.803)
            b.append(3.5098)
            b.append(0.9106)

        # Screening for Nd60
        if Z == 60:
            a.append(0.1665)
            a.append(0.7057)
            a.append(1. - (a[0] + a[1]))
            b.append(23.949)
            b.append(3.5199)
            b.append(0.8486)

        # Screening for Pm61
        if Z == 61:
            a.append(0.1624)
            a.append(0.7133)
            a.append(1. - (a[0] + a[1]))
            b.append(24.598)
            b.append(3.5560)
            b.append(0.8569)

        # Screening for Sm62
        if Z == 62:
            a.append(0.1580)
            a.append(0.7210)
            a.append(1. - (a[0] + a[1]))
            b.append(25.297)
            b.append(3.5963)
            b.append(0.8650)

        # Screening for Eu63
        if Z == 63:
            a.append(0.1538)
            a.append(0.7284)
            a.append(1. - (a[0] + a[1]))
            b.append(26.017)
            b.append(3.6383)
            b.append(0.8731)

        # Screening for Gd64
        if Z == 64:
            a.append(0.1587)
            a.append(0.7024)
            a.append(1. - (a[0] + a[1]))
            b.append(25.497)
            b.append(3.7364)
            b.append(0.9550)

        # Screening for Tb65
        if Z == 65:
            a.append(0.1453)
            a.append(0.7426)
            a.append(1. - (a[0] + a[1]))
            b.append(27.547)
            b.append(3.7288)
            b.append(0.8890)

        # Screening for Dy66
        if Z == 66:
            a.append(0.1413)
            a.append(0.7494)
            a.append(1. - (a[0] + a[1]))
            b.append(28.346)
            b.append(3.7763)
            b.append(0.8969)

        # Screening for Ho67
        if Z == 67:
            a.append(0.1374)
            a.append(0.7558)
            a.append(1. - (a[0] + a[1]))
            b.append(29.160)
            b.append(3.8244)
            b.append(0.9048)

        # Screening for Er68
        if Z == 68:
            a.append(0.1336)
            a.append(0.7619)
            a.append(1. - (a[0] + a[1]))
            b.append(29.990)
            b.append(3.8734)
            b.append(0.9128)

        # Screening for Tm69
        if Z == 69:
            a.append(0.1299)
            a.append(0.7680)
            a.append(1. - (a[0] + a[1]))
            b.append(30.835)
            b.append(3.9233)
            b.append(0.9203)

        # Screening for Yb70
        if Z == 70:
            a.append(0.1267)
            a.append(0.7734)
            a.append(1. - (a[0] + a[1]))
            b.append(31.681)
            b.append(3.9727)
            b.append(0.9288)

        # Screening for Lu71
        if Z == 71:
            a.append(0.1288)
            a.append(0.7528)
            a.append(1. - (a[0] + a[1]))
            b.append(31.353)
            b.append(4.0904)
            b.append(1.0072)

        # Screening for Hf72
        if Z == 72:
            a.append(0.1303)
            a.append(0.7324)
            a.append(1. - (a[0] + a[1]))
            b.append(31.217)
            b.append(4.2049)
            b.append(1.0946)

        # Screening for Ta73
        if Z == 73:
            a.append(0.1384)
            a.append(0.7096)
            a.append(1. - (a[0] + a[1]))
            b.append(30.077)
            b.append(4.2492)
            b.append(1.1697)

        # Screening for W74
        if Z == 74:
            a.append(0.1500)
            a.append(0.6871)
            a.append(1. - (a[0] + a[1]))
            b.append(28.630)
            b.append(4.2426)
            b.append(1.2340)

        # Screening for Re75
        if Z == 75:
            a.append(0.1608)
            a.append(0.6659)
            a.append(1. - (a[0] + a[1]))
            b.append(27.568)
            b.append(4.2341)
            b.append(1.2970)

        # Screening for Os76
        if Z == 76:
            a.append(0.1722)
            a.append(0.6468)
            a.append(1. - (a[0] + a[1]))
            b.append(26.586)
            b.append(4.1999)
            b.append(1.3535)

        # Screening for Ir77
        if Z == 77:
            a.append(0.1834)
            a.append(0.6306)
            a.append(1. - (a[0] + a[1]))
            b.append(25.734)
            b.append(4.1462)
            b.append(1.4037)

        # Screening for Pt78
        if Z == 78:
            a.append(0.2230)
            a.append(0.6176)
            a.append(1. - (a[0] + a[1]))
            b.append(22.994)
            b.append(3.7346)
            b.append(1.4428)

        # Screening for Au79
        if Z == 79:
            a.append(0.2289)
            a.append(0.6114)
            a.append(1. - (a[0] + a[1]))
            b.append(22.864)
            b.append(3.6914)
            b.append(1.4886)

        # Screening for Hg80
        if Z == 80:
            a.append(0.2098)
            a.append(0.6004)
            a.append(1. - (a[0] + a[1]))
            b.append(24.408)
            b.append(3.9643)
            b.append(1.5343)

        # Screening for Tl81
        if Z == 81:
            a.append(0.2708)
            a.append(0.6428)
            a.append(1. - (a[0] + a[1]))
            b.append(20.941)
            b.append(3.2456)
            b.append(1.1121)

        # Screening for Pb82
        if Z == 82:
            a.append(0.2380)
            a.append(0.6308)
            a.append(1. - (a[0] + a[1]))
            b.append(22.987)
            b.append(3.6217)
            b.append(1.2373)

        # Screening for Bi83
        if Z == 83:
            a.append(0.2288)
            a.append(0.6220)
            a.append(1. - (a[0] + a[1]))
            b.append(23.792)
            b.append(3.7796)
            b.append(1.2534)

        # Screening for Po84
        if Z == 84:
            a.append(0.1941)
            a.append(0.6105)
            a.append(1. - (a[0] + a[1]))
            b.append(26.695)
            b.append(4.2582)
            b.append(1.3577)

        # Screening for At85
        if Z == 85:
            a.append(0.1500)
            a.append(0.6031)
            a.append(1. - (a[0] + a[1]))
            b.append(31.840)
            b.append(4.9285)
            b.append(1.4683)

        # Screening for Rn86
        if Z == 86:
            a.append(0.0955)
            a.append(0.6060)
            a.append(1. - (a[0] + a[1]))
            b.append(43.489)
            b.append(5.8520)
            b.append(1.5736)

        # Screening for Fr87
        if Z == 87:
            a.append(0.3192)
            a.append(0.6233)
            a.append(1. - (a[0] + a[1]))
            b.append(20.015)
            b.append(2.9091)
            b.append(0.7207)

        # Screening for Ra88
        if Z == 88:
            a.append(0.2404)
            a.append(0.6567)
            a.append(1. - (a[0] + a[1]))
            b.append(24.501)
            b.append(3.5524)
            b.append(0.8376)

        # Screening for Ac89
        if Z == 89:
            a.append(0.2266)
            a.append(0.6422)
            a.append(1. - (a[0] + a[1]))
            b.append(25.684)
            b.append(3.7922)
            b.append(0.9335)

        # Screening for Th90
        if Z == 90:
            a.append(0.2176)
            a.append(0.6240)
            a.append(1. - (a[0] + a[1]))
            b.append(26.554)
            b.append(4.0044)
            b.append(1.0238)

        # Screening for Pa91
        if Z == 91:
            a.append(0.2413)
            a.append(0.6304)
            a.append(1. - (a[0] + a[1]))
            b.append(25.193)
            b.append(3.6780)
            b.append(0.9699)

        # Screening for U92
        if Z == 92:
            a.append(0.2448)
            a.append(0.6298)
            a.append(1. - (a[0] + a[1]))
            b.append(25.252)
            b.append(3.6397)
            b.append(0.9825)

    # Parameters are in atomic units, thus conversion in natural units
    for i in range(len(b)):
        b[i] = b[i] * ALPHA

    return a, b
