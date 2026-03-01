import math as mt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

#Pedir variables
print("\nDatos Geométricos y del Material")
Sy = float(input("Introduce el Factor de Diseño Von Mises en psi (Sy): "))
D = float(input("Introduce el Diámetro Mayor en mm (D): "))
d = float(input("Introduce el Diámetro Menor en mm (d): "))
r = float(input("Introduce el Radio de Entalle en mm (r): "))
Sut = float(input("Introduce la Resistencia Última (Sut) en psi: "))

print("\nDatos de Carga")
Ma = float(input("Introduce el Momento Alternante (Ma) en psi: "))
Mm = float(input("Introduce el Introduce el Momento Medio (Mm) en psi: "))
Ta = float(input("Introduce el Torque Alternante (Ta) en psi: "))
Tm = float(input("Introduce el Torque Medio (Tm) en psi: "))

print("\nDatos Extras Para las Cuatro Formulas")
Syt = float(input("Resistencia a la Fluentica (Syt) en psi: "))
Se = float(input("Límite de Resistencia a la Fatiga (Se) en psi: "))


#Calcular Kf y Kfs.
#Gráfica del coso Kf, primero el Kt y luego q
X_Kt = np.array([0.02, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.3]) #r/d de Kt (columna)
Curvas_Kt = np.array([1.02, 1.05, 1.10, 1.5, 3]) #Las 5 curvas de Kf (filas)
Y_curvasKt = np.array([
#r/d:0.02, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30
    [2.30, 1.82, 1.65, 1.52, 1.45, 1.39, 1.35, 1.31, 1.28, 1.26, 1.24, 1.22], # curva 1.02 
    [2.58, 1.90, 1.70, 1.58, 1.50, 1.44, 1.39, 1.35, 1.32, 1.29, 1.26, 1.24], # curva 1.05 
    [2.70, 1.98, 1.75, 1.62, 1.55, 1.48, 1.43, 1.39, 1.36, 1.34, 1.32, 1.30], # curva 1.10
    [3, 2.08, 1.84, 1.68, 1.56, 1.48, 1.43, 1.40, 1.37, 1.35, 1.33, 1.3], # curva 1.50
    [3, 2.35, 2.05, 1.85, 1.73, 1.62, 1.55, 1.49, 1.43, 1.39, 1.36, 1.33]  # curva 3.00
])
b_Kt = RegularGridInterpolator( (Curvas_Kt , X_Kt), Y_curvasKt, method='cubic', bounds_error=False, fill_value=None ) #method=cubic->puntos de curvas; bounds_error->predecir el valor

#Grafica de q
X_q = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]) #r,mm (columna)
Curvas_q = np.array([60, 100, 150, 200]) #las 4 curvas (filas)
y_curvasq = np.array([
    #r:0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0
    [0.56, 0.65, 0.70, 0.73, 0.75, 0.77, 0.78, 0.79], # Curva Sut = 60 kpsi
    [0.68, 0.76, 0.80, 0.83, 0.85, 0.87, 0.88, 0.89], # Curva Sut = 100 kpsi
    [0.76, 0.84, 0.88, 0.90, 0.92, 0.93, 0.94, 0.94], # Curva Sut = 150 kpsi
    [0.83, 0.90, 0.93, 0.95, 0.96, 0.97, 0.98, 0.98]  # Curva Sut = 200 kpsi
])
b_q = RegularGridInterpolator( (Curvas_q, X_q), y_curvasq, method='cubic', bounds_error=False, fill_value=None )

#grafica de Kfs
X_Kts = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]) #r/d de Kt (columna)
Curvas_Kts = np.array([1.09, 1.2, 1.33, 2]) #Las 4 curvas de Kf (filas)
Y_curvasKts = np.array([
#r/d:0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
    [1.60, 1.35, 1.21, 1.14, 1.10, 1.08, 1.07], # Curva D/d = 1.09
    [1.95, 1.55, 1.34, 1.23, 1.16, 1.12, 1.10], # Curva D/d = 1.20
    [2.12, 1.65, 1.40, 1.28, 1.20, 1.15, 1.12], # Curva D/d = 1.33
    [2.35, 1.80, 1.48, 1.33, 1.24, 1.18, 1.15]  # Curva D/d = 2.0               
])
b_Kts = RegularGridInterpolator( (Curvas_Kts , X_Kts), Y_curvasKts, method='cubic', bounds_error=False, fill_value=None ) 

#grafica de q_cortante (qc)
X_qc = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]) #r,mm (columna)
Curvas_qc = np.array([60, 100, 150, 200]) #las 4 curvas (filas)
y_curvasqc = np.array([
    #r:0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0
    [0.64, 0.73, 0.78, 0.82, 0.84, 0.86, 0.87, 0.88], # Curva Sut = 60 kpsi
    [0.74, 0.82, 0.86, 0.89, 0.91, 0.92, 0.93, 0.94], # Curva Sut = 100 kpsi
    [0.82, 0.88, 0.91, 0.93, 0.94, 0.95, 0.96, 0.97], # Curva Sut = 150 kpsi
    [0.88, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.98]  # Curva Sut = 200 kpsi
])
b_qc = RegularGridInterpolator( (Curvas_qc, X_qc), y_curvasqc, method='cubic', bounds_error=False, fill_value=None )


#calculos de fatiga:
x_rd = r/d #en mm
c_Dd = D/d #en mm
Sut_kpsi = Sut/1000 #conversión

#calculos de flexión y torsión:
Kt = float( b_Kt([c_Dd, x_rd]) [0]) #flexion
q = float ( b_q([Sut_kpsi, r]) [0])
kf = 1 + q*(Kt-1)

Kts = float( b_Kts([c_Dd, x_rd]) [0]) #torsion
qc = float ( b_qc([Sut_kpsi, r]) [0])
kfs = 1 + qc*(Kts-1)


#Operaciones de las formulas
E_alter = kf * ( (32*Ma) / (mt.pi*pow((d/25.4),3)) )
E_med = kf * ( (32*Mm) / (mt.pi*pow((d/25.4),3)) )
Tau_a = kfs * ( (16*Ta) / (mt.pi*pow((d/25.4),3)) )
Tau_m = kfs * ( (16*Tm) / (mt.pi*pow((d/25.4),3)) )

EPa = pow( ( pow( E_alter, 2) + 3*pow( Tau_a, 2) ) , 0.5) #Esfuerzo prima Alternante
EPm = pow( ( pow( E_med, 2) + 3*pow( Tau_m, 2) ) , 0.5) #Esfuerzo prima Medio
E_max = pow( ( pow((E_med+E_alter) , 2) + 3*pow((Tau_m+Tau_a) , 2)) , 0.5 ) 

Ny = Sy/E_max


#Opreaciones de los Cuatro Factores de Diseños: 
#ln, es 1/n
#goodman
ln_goodman = ( 16/ (mt.pi*pow( (d/25.4), 3)) ) * (    ( (pow( (4*pow( kf*Ma, 2)) + ((3*pow( kfs*Ta, 2))) , 0.5) )/Se )   +    (pow( (4*pow( kf*Mm, 2)) + ((3*pow( kfs*Tm, 2))) , 0.5)/Sut)    )
n_goodman = 1/ln_goodman
#soderberg
ln_soderberg = ( 16/ (mt.pi*pow( (d/25.4), 3)) ) * (    ( (pow( (4*pow( kf*Ma, 2)) + ((3*pow( kfs*Ta, 2))) , 0.5) )/Se )   +    (pow( (4*pow( kf*Mm, 2)) + ((3*pow( kfs*Tm, 2))) , 0.5)/Syt)    )
n_soderberg = 1/ln_soderberg
#gerber
A_gerber = mt.sqrt( (4*pow( kf*Ma, 2)) + (3*pow( kfs*Ta, 2)) )
B_gerber = mt.sqrt( (4*pow( kf*Mm, 2)) + (3*pow( kfs*Tm, 2)) )
ln_gerber = ( (8*A_gerber) / (mt.pi*pow(d/25.4, 3)*Se) )   *   ( 1 + pow( 1 + pow( ( (2*B_gerber*Se) / (A_gerber*Sut )  ) , 2)  , 0.5)  )
n_gerber = 1/ln_gerber
#asme
ln_asme = ( 16/(mt.pi*pow(d/25.4, 3)) ) * pow( (4*pow( (kf*Ma)/Se ,2) ) + (3*pow( (kfs*Ta)/Se ,2) ) + (4*pow( (kf*Mm)/Sy ,2) ) + (3*pow( (kfs*Tm)/Sy ,2) )  , 0.5)
n_asme = 1/ln_asme


#Resultados
print(f"\nCalculos Generales:")
print(f"Sensibilidades a la muesca                                          = {q}, {qc}")
print(f"Factores de concentración de esfuerzo                               = {Kt}, {Kts}")
print(f"Factores de concentración de esfuerzo por fatiga a la flexión (Kf ) = {kf}")
print(f"Factores de concentración de esfuerzo por fatiga a la torsión (Kfs) = {kfs}")

print(f"\nEsfuerzo Alternante (σa)         = {E_alter}")
print(f"Esfuerzo Medio (σm)                = {E_med}")
print(f"Esfuerzo Cortante Alternante (τa)  = {Tau_a}")
print(f"Esfuerzo Cortante Medio (τm)       = {Tau_m}")

print(f"\nEsfuerzo Prima Alternandes de Von Mises (σ'a) = {EPa}")
print(f"Esfuerzo Prima Medio de Von Mises (σ'm)         = {EPm}")
print(f"Esfuerzo Máximo de Von Mises (σ'max)            = {E_max}")
print(f"\nLa Resistencia a la Influencia (Ny)           = {Ny}")

print(f"\n\nLos Cuatro Factores de Diseños:")
print(f"\nCríterio de ED-Goodman:")
print(f"1/n = {ln_goodman}")
print(f"n = {n_goodman}")
print(f"\nCríterio de ED-Soderberg:")
print(f"1/n = {ln_soderberg}")
print(f"n = {n_soderberg}")
print(f"\nCríterio de ED-Gerber:")
print(f"1/n = {ln_gerber}")
print(f"n = {n_gerber}")
print(f"\nCríterio de ED-ASME elíptica:")
print(f"1/n = {ln_asme}")
print(f"n = {n_asme}")



#graficar las graficas daaaaaa
#axs, guarda los 4 recuadros;   plt.subplots, ventana emergente;   figsize, tamaño en in;   suptitle, el titulo centrado
fig, axs = plt.subplots(2, 2, figsize=(14, 10) )
fig.suptitle( 'Demostración visual de los Factores', fontsize=16, fontweight='bold' )

#lista de 100 numeros mediante un intervalo, para que se vea bonito, po'quesi
#for Dd in Curvas_Kt, repite el proceso de las 5 curvas;   ly_kt, lista de compresión, toma el valor de la curva y calcula para los 100
xl_kt = np.linspace(0.02, 0.3, 100) #grafica Kt
for Dd in Curvas_Kt: #las 5 curvas de kt
    yl_kt = [b_Kt([Dd, x])[0] for x in xl_kt ]
    axs[0, 0].plot(xl_kt, yl_kt, label=f'D/d = {Dd}', linewidth=2 )

#punto rojo en la intersección
#marker='o', en lugar de lineas, solo el punto;  markersize, tamaño del punto;  zorder, capas 
axs[0, 0].plot(x_rd, Kt, marker='o', color='red', markersize=5, zorder=5, label='Valor Kt')

#estetico
#.set, detos de la grafica;   .grid, enciende la cuadricula;   linestyle, hace que sea punteada;   alpha, lo hace transparente
#.lengend, muestra el cuadrito de información (color-D/d)
axs[0, 0].set_title('Kt - Flexión')
axs[0, 0].set_xlabel('r/d')
axs[0, 0].set_ylabel('Kt')
axs[0, 0].grid(True, linestyle='--', alpha=0.7 )
axs[0, 0].legend()


#ahora la grafica q
xl_q = np.linspace(0.5, 4, 100) #si flota el punto, es porque r es <0.5mm
for Sut in Curvas_q: #las 4 curvas de q
    yl_q = [b_q([Sut, x])[0] for x in xl_q ]
    axs[0, 1].plot(xl_q, yl_q, label=f'Sut = {Sut} kpsi', linewidth=2 )

axs[0, 1].plot(r, q, marker='o', color='red', markersize=5, zorder=5, label='Valor q') #el coso rojo

#estetico
axs[0, 1].set_title('q - Flexión')
axs[0, 1].set_xlabel('r (mm)')
axs[0, 1].set_ylabel('q')
axs[0, 1].grid(True, linestyle='--', alpha=0.7 )
axs[0, 1].legend()


#grafica Kts
xl_kts = np.linspace(0.02, 0.3, 100)
for Dd in Curvas_Kts: #las 5 curvas de kt
    yl_kts = [b_Kts([Dd, x])[0] for x in xl_kts ]
    axs[1, 0].plot(xl_kts, yl_kts, label=f'D/d = {Dd}', linewidth=2 )

axs[1, 0].plot(x_rd, Kts, marker='o', color='red', markersize=5, zorder=5, label='Valor Kts')

#estetico
axs[1, 0].set_title('Kts - Torsión')
axs[1, 0].set_xlabel('r/d')
axs[1, 0].set_ylabel('Kts')
axs[1, 0].grid(True, linestyle='--', alpha=0.7 )
axs[1, 0].legend()


#ahora la grafica qcortante
xl_qc = np.linspace(0.5, 4, 100)
for Sut in Curvas_qc: #las 4 curvas de q
    yl_qc = [b_qc([Sut, x])[0] for x in xl_qc ]
    axs[1, 1].plot(xl_qc, yl_qc, label=f'Sut = {Sut} kpsi', linewidth=2 )

axs[1, 1].plot(r, qc, marker='o', color='red', markersize=5, zorder=5, label='Valor qc')

#estetico
axs[1, 1].set_title('q cortante - Torsión')
axs[1, 1].set_xlabel('r (mm)')
axs[1, 1].set_ylabel('q cortante')
axs[1, 1].grid(True, linestyle='--', alpha=0.7 )
axs[1, 1].legend()


#limpieza visual
#plt.tight_layout, automatico, agusta los margenes para que no se encimen,   plt.show, comando de ejecutión, detiene el código y espera para que se vean asi bonito jeje
plt.tight_layout()
plt.show()

