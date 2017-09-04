# -*- coding: utf-8 -*-
"""
Редактор Spyder
 
Это временный скриптовый файл.
"""



# Загрузка библиотек необходимых для отрисовки графиков
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d 
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy.optimize import broyden1
from scipy.optimize import fsolve
from scipy.optimize import anderson
#%matplotlib inline 

""" 
# Исходные данные для проведения расчета
v_s_l = 2          # приведенная скорость жидкости, м3/с
v_s_g = 5            # приведенная скорость газа, м3/с
P_pipe = 40            # давление в трубе, МПа
T_pipe = 350          # температура на устье скважины, К
L_pipe = 2000          # длина трубы, м
d_pipe = 0.8     # диаметр трубы, м
rho_l_PT = 1000      # плотность жидкости при термобарических условиях, кг/м3
rho_g_PT = 10      # плотность газа при термобарических условиях, кг/м3
mu_l_PT = 0.001    # вязкость жидкости при термобарических условиях, Па*с
mu_g_PT = 0.00002    # вязкость газа при термобарических условиях, Па*с
angle = 0     # угол наклона трубы, градус
rough = 0.0001   # абсолютная шероховатость трубы, м
gravity = 9.81  #ускорение свободного падения, м/с2
sigma_l = 0.03 #поверхностное натяжение на границе жидкости и газа, Н/м

""" 

class grad_Xiao:     
    # Расчет предварительных величин
    def __init__(self):
        """
        конструктор класса, задает исходные значения 
        """
        # Исходные данные для проведения расчета
        self.v_s_l = 2          # приведенная скорость жидкости, м3/с
        self.v_s_g = 5            # приведенная скорость газа, м3/с
        self.P_pipe = 40            # давление в трубе, МПа
        self.T_pipe = 350          # температура на устье скважины, К
        self.L_pipe = 2000          # длина трубы, м
        self.d_pipe = 0.8     # диаметр трубы, м
        self.rho_l_PT = 1000      # плотность жидкости при термобарических условиях, кг/м3
        self.rho_g_PT = 10      # плотность газа при термобарических условиях, кг/м3
        self.mu_l_PT = 0.001    # вязкость жидкости при термобарических условиях, Па*с
        self.mu_g_PT = 0.00002    # вязкость газа при термобарических условиях, Па*с
        self.angle = 0     # угол наклона трубы, градус
        self.rough = 0.0001   # абсолютная шероховатость трубы, м
        self.gravity = 9.81  #ускорение свободного падения, м/с2
        self.sigma_l = 0.03 #поверхностное натяжение на границе жидкости и газа, Н/м   
        
        self.Area_pipe_m =  math.pi*self.d_pipe**2/4  # диаметр трубы
    
    def calc_Area_pipe(self, d_pipe):      # расчет поперечной площади трубы, м2
        return math.pi*d_pipe**2/4
     
    def calc_dim_rough(self, rough,d_pipe):       # расчет относительной шероховатости трубы
        return rough/1000/d_pipe
    # Расчет скоростей
     
     
    # Расчет приведенных скоростей - Calculate Superficial Velocity
    def calc_lyambda_l_PT(self, v_s_l, v_s_g):      # расчет рсходного объемного содержания жидкости при термобарических условиях, м3/сут
        return v_s_l/(v_s_l+v_s_g)
    
     
    # Расчет скорости смеси - Calculate Superficial Velocity
    def calc_v_s_mix(self, v_s_l, v_s_g):
        return v_s_g+v_s_l
     
    # Расчет действительных скоростей - Calculate Actual Velocity
    def calc_v_a_g(self, v_s_l, v_s_g):    # Actual Gas Velocity
        lyambda_l_PT=v_s_l/(v_s_l+v_s_g)
        return v_s_g/(1-lyambda_l_PT)
     
    def calc_v_a_l(self, v_s_l, v_s_g):    # Actual Liquid Velocity
        lyambda_l_PT=v_s_l/(v_s_l+v_s_g)
        return v_s_l/lyambda_l_PT
       
    # Расчет скорости проскальзывания газа
    def calc_v_slip(self,v_s_l, v_s_g):
        v_a_g=self.calc_v_a_g(v_s_l, v_s_g)
        v_a_l=self.calc_v_a_l(v_s_l, v_s_g)
        return v_a_g-v_a_l
     
    # Расчет градиента давления - Pressure Gradient Prediction
     
    # Расчет чисел Рейнольдса для жидкости и газа
    def calc_N_Re_l(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):     # Расчет числа Рейнольдса для жидкости
        return rho_l_PT*d_pipe*v_s_l/mu_l_PT
     
    def calc_N_Re_g(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT):     # Расчет числа Рейнольдса для газа
        return rho_g_PT*d_pipe*v_s_g/mu_g_PT
     
    # Расчет гидравлического коэффициента потерь на трение - Calculate Friction Factor
    # Зарубежная корреляция для расчета коэффициента гиравлических потерь на трение - Moody Chart
    # Принимаем допущение, что при Re>2000 - режим турбулентный, нет переходного режима, так как не определены коэффициенты
    # для friction factor in transition region (friction factor=C*Re^n)
    def calc_f_factor_l(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):     # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.calc_N_Re_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT) < 2000:
        # режим ламинарный
            return 64*self.calc_N_Re_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT)**(-0.2)
     
    def calc_f_factor_g(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT):     # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.calc_N_Re_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)<2000:
        # режим ламинарный
            return 64*self.calc_N_Re_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)**(-0.2)
           
            
    def calc_pres_drop_s_l(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):        # Оценка градиента давления для жидкости, для приведенной скорости
        f_factor_l=self.calc_f_factor_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
        return -f_factor_l*rho_l_PT*v_s_l**2/(2*d_pipe)
       
    def calc_pres_drop_s_g(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT):        # Оценка градиента давления для газа, для приведенной скорости
        f_factor_g=self.calc_f_factor_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
        return -f_factor_g*rho_g_PT*v_s_g**2/(2*d_pipe)
           
    # Расчет параметров, определяющих границы перехода между режимами течения
     
    # Расчет вспомогательных параметров для определения безрамерной высоты жидкости - Calculate Dimensionless Liquid Height
    
    def calc_dim_h_l_fact(self,dim_h_l):
        if dim_h_l>=1:
            return 0.9999
        elif dim_h_l<=0:
            return 0.0001
        else:
            return dim_h_l
    
    def calc_dim_A_l(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        b=2*dim_h_l-1
        a=math.pi - math.acos(b)
        c=1-b*b
        return 0.25*(a+b*c**0.5)
     
    def calc_dim_A_g(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        b=2*dim_h_l-1    
        a=math.acos(b)
        c=1-b**2
        return 0.25*(a-b*c**0.5)
     
    def calc_dim_S_l(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        return math.pi - math.acos(2*dim_h_l-1)
     
    def calc_dim_S_i(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        return (1-(2*dim_h_l-1)**2)**0.5
       
    def calc_dim_S_g(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        return math.acos(2*dim_h_l-1)
       
    def calc_dim_v_l(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        return (dim_A_g/dim_A_l+1)
       
    def calc_dim_v_g(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        return 1+dim_A_l/dim_A_g
       
    def calc_dim_d_l(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        dim_S_l=math.pi - math.acos(2*dim_h_l-1)
        return 4*dim_A_l/dim_S_l
       
    def calc_dim_d_g(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
        dim_S_g=math.acos(2*dim_h_l-1)
        dim_S_i=(1-(2*dim_h_l-1)**2)**0.5
        return 4*dim_A_g/(dim_S_i+dim_S_g)
       
    def calc_n_power_coeff(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):
        N_Re_l=self.calc_N_Re_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
        if N_Re_l<2000:
            return 1
        else:
            return 0.2
     
    def calc_m_power_coeff(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT):
        N_Re_g=self.calc_N_Re_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)    
        if N_Re_g<2000:
            return 1
        else:
            return 0.2
       
    def calc_x_large(self,rho_l_PT,rho_g_PT,v_s_l,v_s_g,d_pipe,mu_l_PT,mu_g_PT):      # ratio of the superficial liquid and superficial gas frictional pressure gradients
        pres_drop_s_l=self.calc_pres_drop_s_l(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
        pres_drop_s_g=self.calc_pres_drop_s_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
        return ((-pres_drop_s_l)/(-pres_drop_s_g))**0.5
     
    def calc_y_large(self,rho_l_PT,rho_g_PT,v_s_g,d_pipe,mu_g_PT,gravity,angle):      # to introduce effect of inclination angle on the height of the stratified liquid layer
        pres_drop_s_g=self.calc_pres_drop_s_g(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
        return (rho_l_PT-rho_g_PT)*gravity*math.sin(math.pi*angle/180)/(-pres_drop_s_g)
     
    def calc_uniq_func_dim_height(self,dim_h_l):
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        a1=self.calc_x_large(self.rho_l_PT,self.rho_g_PT,self.v_s_l,self.v_s_g,self.d_pipe,self.mu_l_PT,self.mu_g_PT)
        a2=self.calc_y_large(self.rho_l_PT,self.rho_g_PT,self.v_s_g,self.d_pipe,self.mu_g_PT,self.gravity,self.angle)
        b1=self.calc_dim_A_l(dim_h_l)
        b2=self.calc_dim_A_g(dim_h_l)
        c1=self.calc_dim_S_l(dim_h_l)
        c2=self.calc_dim_S_i(dim_h_l)
        c3=self.calc_dim_S_g(dim_h_l)
        d1=self.calc_dim_v_l(dim_h_l)
        d2=self.calc_dim_v_g(dim_h_l)
        e1=self.calc_dim_d_l(dim_h_l)
        e2=self.calc_dim_d_g(dim_h_l)
        f1=self.calc_n_power_coeff(self.rho_l_PT,self.v_s_l,self.d_pipe,self.mu_l_PT)
        f2=self.calc_m_power_coeff(self.rho_g_PT,self.v_s_g,self.d_pipe,self.mu_g_PT)
        a11=a1*a1
        b11=d1*e1
        c11=b11**(-f1)
        d11=d1*d1
        e11=c1/b1
        f11=d2*e2
        g11=f11**(-f2)
        h11=d2*d2
        i11=c3/b2+c2/b1+c2/b2
        return  a11*c11*d11*e11-g11*h11*i11+4*a2
     
     
    #определение струтуры потока
     
    def calc_f_large(self,gravity,angle,d_pipe,rho_l_PT,rho_g_PT,v_s_l,v_s_g):
        return (rho_g_PT/(rho_l_PT-rho_g_PT))**0.5*v_s_g/(d_pipe*gravity*math.cos(math.pi*angle/180))**0.5
     
    def calc_k_large(self):
        N_Re_l=self.calc_N_Re_l(self.rho_l_PT,self.v_s_l,self.d_pipe,self.mu_l_PT)
        f_large=self.calc_f_large(self.gravity,self.angle,self.d_pipe,self.rho_l_PT,self.rho_g_PT,self.v_s_l,self.v_s_g)
        return N_Re_l**0.5*f_large
     
    def calc_t_large(self):
        pres_drop_s_l=self.calc_pres_drop_s_l(self.rho_l_PT,self.v_s_l,self.d_pipe,self.mu_l_PT)
        return (-pres_drop_s_l/((self.rho_l_PT-self.rho_g_PT)*self.gravity*math.cos(math.pi*self.angle/180)))**0.5
    
    def calc_trans_A(self,dim_h_l):             #Stratified to Nonstratified Transition criterea (если критерий больше либо равен единице, то Nonstratified)
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_v_g=self.calc_dim_v_g(dim_h_l)
        dim_S_i=self.calc_dim_S_i(dim_h_l)
        dim_A_g=self.calc_dim_A_g(dim_h_l)
        f_large=self.calc_f_large(self.gravity,self.angle,self.d_pipe,self.rho_l_PT,self.rho_g_PT,self.v_s_l,self.v_s_g)
        return f_large**2*(1/(1-dim_h_l)**2*(dim_v_g**2*dim_S_i/dim_A_g))
     
    def calc_trans_B(self,dim_h_l):        #Intermittent or Dispersed_bubble to Annular Transition criterea (если критерий меньше либо равен 0,35, то Annular)
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        return dim_h_l
     
    def calc_trans_C(self,dim_h_l):         #Stratified-Smooth to Stratified-Wavy Transition criterea (если больше либо равно нулю, то Stratified-Smooth)
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_v_g=self.calc_dim_v_g(dim_h_l)
        dim_v_l=self.calc_dim_v_l(dim_h_l)
        s=0.01
        k_large=self.calc_k_large()
        return k_large**2-2/(dim_v_l**0.5*dim_v_g*s**0.5)
     
    def calc_trans_D(self,dim_h_l):         #Intermittent to Dispersed-Bubble Transition criterea (если критерий больше либо равен нулю, то Dispersed-Bubble)
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        dim_v_l=self.calc_dim_v_l(dim_h_l)
        dim_S_i=self.calc_dim_S_i(dim_h_l)
        dim_A_g=self.calc_dim_A_g(dim_h_l)
        dim_d_l=self.calc_dim_d_l(dim_h_l)
        n_power_coeff=self.calc_n_power_coeff(self.rho_l_PT,self.v_s_l,self.d_pipe,self.mu_l_PT)
        t_large=self.calc_t_large()
        return t_large**2-(8*dim_A_g/(dim_S_i*dim_v_l**2*(dim_v_l*dim_d_l)**(-n_power_coeff)))
       
    def calc_stratified_regime(self,dim_h_l):   #если режим определен, что stratified, то выбираем wavy or smooth
                                    #1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        trans_C=self.calc_trans_C(dim_h_l)
        if trans_C<0:
            return 3
        else:
            return 2
       
    def calc_bubble_slug_regime(self,dim_h_l):  #check dispersed-bubble or intermittent(i.e. slug or plug)
                                    #1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        trans_D=self.calc_trans_D(dim_h_l)
        if trans_D<0:
            return 5
        else:
            return 4
     
    def calc_nonstratified_regime(self,dim_h_l):#если режим определен, что nonstratified, то выбираем между annular or intermittent or dispersed-bubble
                                    #1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        trans_A=self.calc_trans_A(dim_h_l)
        trans_B=self.calc_trans_B(dim_h_l)
        bubble_slug_regime=self.calc_bubble_slug_regime(dim_h_l)
        if trans_A>=1 and trans_B<=0.35:
    #        return 1        #Annular Regime будет наблюдаться в трубах малого диаметра, где есть эффект капилляра
            return 2        #Мы принимаем, что у нас Stratified-Wavy Regime, в связи с тем что площадь стенки
        #очень большая жидкость из-за гравитации будет стекать, не наблюдаетяс эффект капилляра в большой трубе
        else:
            return bubble_slug_regime     
     
    def calc_flow_structure(self,dim_h_l): #1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        trans_A=self.calc_trans_A(dim_h_l)
        nonstratified_regime=self.calc_nonstratified_regime(dim_h_l)
        stratified_regime=self.calc_stratified_regime(dim_h_l)
        if trans_A<1:
            return stratified_regime
        else:
            return nonstratified_regime
     
    

     
    #for xx in t :
    #    yy.append(calc_uniq_func_dim_height(xx))
    #plt.plot(t,yy)   
     
    #print("dim_h_l = ",dim_h_l)
        
    #print(calc_v_s_g(Q_g_PT,d_pipe))
    #print(calc_v_s_l(Q_o_PT,Q_w_PT,d_pipe))
    #print(calc_v_mix(Q_o_PT,Q_w_PT,Q_g_PT,d_pipe))
     
    #print(calc_lyambda_l_PT(Q_o_PT,Q_w_PT,Q_g_PT))
     
    #print(calc_f_large(rho_g_PT,gravity,angle,d_pipe,rho_o_PT,rho_w_PT,Q_w_PT,Q_o_PT))
    #print(calc_k_large())
    #print (calc_t_large())
    #print(calc_trans_A())
    #print(calc_trans_B())
    #print(calc_trans_C())
    #print(calc_trans_D())
    
    
    
    
    
    #Liquid Holdup Prediction
    
    #Stratified flow
    def calc_H_L_strat(self,dim_h_l):
        h_l=self.calc_dim_h_l_fact(dim_h_l)
        return 0.5 + (4/math.pi)*(h_l - 0.5)*(h_l-h_l**2)**0.5 + (1/math.pi)*math.asin(2*h_l-1)
    
    #Annular flow
    def calc_N_We_s_g_an(self,rho_g_PT, v_s_g, d_pipe, sigma_l):   #Число Вебера в ядре газа
        return rho_g_PT*v_s_g**2*d_pipe/sigma_l
    
    def calc_N_Re_s_g_an(self,rho_g_PT, v_s_g, d_pipe, mu_g_PT):     #Число Рейнольдса для газа в ядре газа
        return rho_g_PT*v_s_g*d_pipe/mu_g_PT
    
    def calc_N_Re_s_l_an(self,rho_l_PT, v_s_l, d_pipe, mu_l_PT):     #Число Рейнольдса для жидкости в ядре газа
        return rho_l_PT*v_s_l*d_pipe/mu_l_PT
    
    def calc_N_Fr_s_g_an(self,v_s_g, d_pipe):
        return v_s_g/(self.gravity*d_pipe)**0.5
    
    def calc_f_E_an(self,v_s_g, v_s_l):    #Доля увлеченной жидкости в газовом ядре
        N_We_s_g_an = self.calc_N_We_s_g_an(self.rho_g_PT, v_s_g, self.d_pipe, self.sigma_l)
        N_Re_s_g_an = self.calc_N_Re_s_g_an(self.rho_g_PT, self.v_s_g, self.d_pipe, self.mu_g_PT)
        N_Re_s_l_an = self.calc_N_Re_s_l_an(self.rho_l_PT, self.v_s_l, self.d_pipe, self.mu_l_PT)
        N_Fr_s_g_an = self.calc_N_Fr_s_g_an(v_s_g, self.d_pipe)
        c1 = N_We_s_g_an**(1.8)
        c2 = N_Fr_s_g_an**(-0.92)
        c3 = N_Re_s_l_an**0.7
        c4 = N_Re_s_g_an**(-1.24)
        c5 = (self.rho_l_PT/self.rho_g_PT)
        c55 = c5**0.38
        c6 = (self.mu_l_PT/self.mu_g_PT)
        c66 = c6**0.97
        c = 0.003*c1*c2*c3*c4*c55*c66
        return c/(c+1)
    
    def calc_H_L_gas_core(self,v_s_g, v_s_l):    #Gas Core Liquid Holdup
        f_E = self.calc_f_E_an(v_s_g, v_s_l)
        return v_s_l*f_E/(v_s_g+v_s_l*f_E)
    
    def calc_v_s_g_modified(self,v_s_g, v_s_l):  #модифицированные параметры для ядра, исходя из предпложения, что поток в яжре гомогенный
        f_E = self.calc_f_E_an(v_s_g, v_s_l)
        return v_s_g + f_E*v_s_l
    
    def calc_v_s_l_modified(self,v_s_g, v_s_l):  
        f_E = self.calc_f_E_an(v_s_g, v_s_l)
        return v_s_l - f_E*v_s_l
    
    def calc_rho_g_modified(self,v_s_g, v_s_l):
        f_E = self.calc_f_E_an(v_s_g, v_s_l)
        v_s_g_modified = self.calc_v_s_g_modified(v_s_g, v_s_l)
        return (self.rho_g_PT*v_s_g + self.rho_l_PT*f_E*v_s_l)/v_s_g_modified
    
    def calc_mu_g_modified(self,v_s_g, v_s_l):
        f_E = self.calc_f_E_an(v_s_g, v_s_l)
        v_s_g_modified = self.calc_v_s_g_modified(v_s_g, v_s_l)
        return (self.mu_g_PT*v_s_g + self.mu_l_PT*f_E*v_s_l)/v_s_g_modified
    
    def calc_N_Re_l_modified(self,rho_l_PT,v_s_l,v_s_g,d_pipe,mu_l_PT):     # Расчет числа Рейнольдса для жидкости с учетом модифицированных параметров
        v_s_l_modified = self.calc_v_s_l_modified(v_s_g, v_s_l)
        return rho_l_PT*d_pipe*v_s_l_modified/mu_l_PT
     
    def calc_N_Re_g_modified(self,rho_g_PT,v_s_g,v_s_l,d_pipe,mu_g_PT):     # Расчет числа Рейнольдса для газа с учетом модифицированных параметров
        v_s_g_modified = self.calc_v_s_g_modified(v_s_g, v_s_l)
        rho_g_PT_modified = self.calc_rho_g_modified(v_s_g, v_s_l)
        mu_g_PT_modified = self.calc_mu_g_modified(v_s_g, v_s_l)
        return rho_g_PT_modified*d_pipe*v_s_g_modified/mu_g_PT_modified
    
    def calc_f_factor_l_modified(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):     # Расчет гидравлического коэффициента потерь на трение для жидкости с учетом модифицированных параметров
        if self.calc_N_Re_l_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT) < 2000:
        # режим ламинарный
            return 64*self.calc_N_Re_l_modified(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_l_modified(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT)**(-0.2)
     
    def calc_f_factor_g_modified(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT):     # Расчет гидравлического коэффициента потерь на трение для жидкости с учетом модифицированных параметров
        if self.calc_N_Re_g_modified(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT)<2000:
        # режим ламинарный
            return 64*self.calc_N_Re_g_modified(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_g_modified(self,rho_g_PT,v_s_g,d_pipe,mu_g_PT)**(-0.2)
    
    
    def calc_pres_drop_s_l_modified(self,rho_l_PT,v_s_l,d_pipe,mu_l_PT):        # Оценка градиента давления для жидкости, для приведенной скорости с учетом модифицированных параметров
        f_factor_l_modified=self.calc_f_factor_l_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
        v_s_l_modified = v_s_l
    #    v_s_l_modified = calc_v_s_l_modified(v_s_g, v_s_l)
        return -f_factor_l_modified*rho_l_PT*v_s_l_modified**2/(2*d_pipe)
       
    def calc_pres_drop_s_g_modified(self,rho_g_PT,v_s_g,v_s_l,d_pipe,mu_g_PT):        # Оценка градиента давления для газа, для приведенной скорости с учетом модифицированных параметров
        v_s_g_modified = self.calc_v_s_g_modified(v_s_g, v_s_l)
        rho_g_PT_modified = self.calc_rho_g_modified(v_s_g, v_s_l)
#        mu_g_PT_modified = self.calc_mu_g_modified(v_s_g, v_s_l)
        f_factor_g_modified=self.calc_f_factor_g_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
        return -f_factor_g_modified*rho_g_PT_modified*v_s_g_modified**2/(2*d_pipe)
    
    # С этого момента расчет холдапа в пленке в Annular режиме, как писал выше опускаем данный пункт, либо принимаем холдап равный холдапу в ядре
    
    #def calc_H_L_f_an_fact_modified(H_L_f_an):
    #    if H_L_f_an>=1:
    #        return 0.001
    #    elif H_L_f_an<=0.00001:
    #        return 0.001
    #    else:
    #        return H_L_f_an
    
    #def calc_dim_A_l_modified(H_L_f_an):  #!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    b=2*H_L_f_an-1
    #    a=math.pi - math.acos(b)
    #    c=1-b*b
    #    return 0.25*(a+b*c**0.5)
     
    #def calc_dim_A_g_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    b=2*dim_h_l-1    
    #    a=math.acos(b)
    #    c=1-b**2
    #    return 0.25*(a-b*c**0.5)
     
    #def calc_dim_S_l_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    return math.pi - math.acos(2*dim_h_l-1)
     
    #def calc_dim_S_i_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    return (1-(2*dim_h_l-1)**2)**0.5
       
    #def calc_dim_S_g_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    return math.acos(2*dim_h_l-1)
       
    #def calc_dim_v_l_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    return (dim_A_g/dim_A_l+1)
       
    #def calc_dim_v_g_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    return 1+dim_A_l/dim_A_g
       
    #def calc_dim_d_l_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    dim_S_l=math.pi - math.acos(2*dim_h_l-1)
    #    return 4*dim_A_l/dim_S_l
       
    #def calc_dim_d_g_modified(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
    #    dim_S_g=math.acos(2*dim_h_l-1)
    #    dim_S_i=(1-(2*dim_h_l-1)**2)**0.5
    #    return 4*dim_A_g/(dim_S_i+dim_S_g)
    
    #def calc_n_power_coeff_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT):
    #    N_Re_l_modified=calc_N_Re_l_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
    #    if N_Re_l_modified<2000:
    #        return 1
    #    else:
    #        return 0.2
     
    #def calc_m_power_coeff_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT):
    #    N_Re_g_modified=calc_N_Re_g_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT)    
    #    if N_Re_g_modified<2000:
    #        return 1
    #    else:
    #        return 0.2
    
    #def calc_x_large_modified(rho_l_PT,v_s_l,v_s_g,d_pipe,mu_l_PT):      # ratio of the superficial liquid and superficial gas frictional pressure gradients
    #    pres_drop_s_l=calc_pres_drop_s_l_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
    #    pres_drop_s_g=calc_pres_drop_s_g_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
    #    return ((-pres_drop_s_l)/(-pres_drop_s_g))**0.5
     
    #def calc_y_large_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT,gravity,angle):      # to introduce effect of inclination angle on the height of the stratified liquid layer
    #    pres_drop_s_g=calc_pres_drop_s_g_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
    #    rho_g_PT_modified = calc_rho_g_modified(v_s_g, v_s_l)
    #    return (rho_l_PT-rho_g_PT_modified)*gravity*math.sin(math.pi*angle/180)/(-pres_drop_s_g)
     
    #def calc_H_L_f_an(H_L_f_an):
    #    H_L_f_an=calc_H_L_f_an_fact_modified(H_L_f_an)
    #    a1=calc_x_large_modified(rho_l_PT,v_s_l,v_s_g,d_pipe,mu_l_PT)
    #    a2=calc_y_large_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT,gravity,angle)
    #    b1=calc_dim_A_l_modified(H_L_f_an)
    #    b2=calc_dim_A_g_modified(H_L_f_an)
    #    c1=calc_dim_S_l_modified(H_L_f_an)
    #    c2=calc_dim_S_i_modified(H_L_f_an)
    #    c3=calc_dim_S_g_modified(H_L_f_an)
    #    d1=calc_dim_v_l_modified(H_L_f_an)
    #    d2=calc_dim_v_g_modified(H_L_f_an)
    #    e1=calc_dim_d_l_modified(H_L_f_an)
    #    e2=calc_dim_d_g_modified(H_L_f_an)
    #    f1=calc_n_power_coeff_modified(rho_l_PT,v_s_l,d_pipe,mu_l_PT)
    #    f2=calc_m_power_coeff_modified(rho_g_PT,v_s_g,d_pipe,mu_g_PT)
    #    a11=a1*a1
    #    b11=d1*e1
    #    c11=b11**(-f1)
    #    d11=d1*d1
    #    e11=c1/b1
    #    f11=d2*e2
    #    g11=f11**(-f2)
    #    h11=d2*d2
    #    i11=c3/b2+c2/b1+c2/b2
    #    return  a11*c11*d11*e11-g11*h11*i11+4*a2
    
    #from scipy.optimize import 
    
    #H_L_f_an = fsolve(calc_H_L_f_an, 0.0001) #Liquid Holdup in Annular Film
    
    #def calc_H_L_total_an(v_s_g, v_s_l):    #Total Annular Flow Liquid Holdup
    #    x = calc_H_L_f_an_fact_modified(H_L_f_an)
    #    H_L_gas_core = calc_H_L_gas_core(v_s_g, v_s_l)
    #    return x + (1-x)*H_L_gas_core
    
    #Dispersed Bubble Flow
    def calc_H_L_dispers(self,v_s_g, v_s_l):
        return v_s_l/(v_s_l+v_s_g)
    
    #Intermittent Flow
    def calc_v_s(self,v_s_g, v_s_l):         #Slug Velosity - средняя скорость гомогенной смеси в Slug Region (fig 4.14)
        return self.calc_v_s_mix(v_s_l, v_s_g)
    
    def calc_Re_mix_slug(self,v_s_g, v_s_l, rho_l_PT, mu_l_PT, d_pipe):      #расчет числа Рейнольдса для смеси при пробковом режиме, согласно  Fabre(1994).27961-MS.Advancements in Two-Phase Slug Flow Modeling 
        v_s = self.calc_v_s(v_s_g, v_s_l)
        return v_s*rho_l_PT*d_pipe/mu_l_PT
    
    def calc_c_o(self,v_s_g, v_s_l):         #Коэффициент, учитывающий насколько максимальная скорость больше средней для фронта Slug, согласно  Fabre(1994).27961-MS.Advancements in Two-Phase Slug Flow Modeling 
        Re_mix_slug_crit = 1000
        Re_mix_slug = self.calc_Re_mix_slug(v_s_g, v_s_l, self.rho_l_PT, self.mu_l_PT, self.d_pipe)
        return 2.27/(1+(Re_mix_slug/Re_mix_slug_crit)**2)+1.2/(1+(Re_mix_slug_crit/Re_mix_slug)**2)
    
    def calc_v_d_TB(self,v_s_g, v_s_l):      #Taylor Bubble Drift Velocity
        a = (self.gravity*self.d_pipe)**0.5
        return 0.54*a*math.cos(math.pi*self.angle/180)+0.35*a**0.5*math.sin(math.pi*self.angle/180)
    
    def calc_v_T_B(self,v_s_g, v_s_l):       #Translational Veosity на границе slug региона и региона пуырька Тейлора
        c_o = self.calc_c_o(v_s_g, v_s_l)
        v_d_TB = self.calc_v_d_TB(v_s_g, v_s_l)
        v_s = self.calc_v_s(v_s_g, v_s_l)
        return c_o*v_s + v_d_TB
    
    def calc_v_g_L_S(self,v_s_g, v_s_l):         #Gas Velocity in Liquid Slug
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        v_s = self.calc_v_s(v_s_g, v_s_l)
        c = (self.gravity*self.sigma_l*(self.rho_l_PT-self.rho_g_PT)/self.rho_l_PT**2)**0.25
        return v_s + 1.53*c*H_L_L_S**0.1*math.sin(math.pi*self.angle/180)
        
    def calc_v_L_L_S(self,v_s_g, v_s_l):         #Liquid Velocity in Liquid Slug
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        v_g_L_S = self.calc_v_g_L_S(v_s_g, v_s_l)
        return (v_s_l + v_s_g - v_g_L_S*(1 - H_L_L_S))/H_L_L_S
    
    def calc_H_L_L_S(self,v_s_g, v_s_l):         #Slug Liquid Holdup in Slug Region (fig 4.14)
        v_s = self.calc_v_s(v_s_g, v_s_l)
        c = (v_s/8.66)**1.39
        return 1/(1+c)
    
    def calc_H_L_S_U(self,v_s_g, v_s_l):         #Average Liquid Holdup
        v_T_B = self.calc_v_T_B(v_s_g, v_s_l)
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        v_g_L_S = self.calc_v_g_L_S(v_s_g, v_s_l)
        return (v_T_B*H_L_L_S + v_g_L_S*(1-H_L_L_S) - v_s_g)/v_T_B
    
    #Определение Liquid Holdup in Taylor Bubble
    def calc_H_L_T_B(self,h_f_d):        #Вычисление холдапа жидкоти в Тейлор бабл регионе из предположения размера пленки, как для стратифайд режима
        return 0.5 + (4/math.pi)*(h_f_d - 0.5)*(h_f_d-h_f_d**2)**0.5 + (1/math.pi)*math.asin(2*h_f_d-1)
    
    def calc_v_L_T_B(self,v_s_g, v_s_l):         #Actual film velocity in Taylor Bubble Region
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        H_L_T_B = self.calc_H_L_T_B(self.h_f_d)
        v_T_B = self.calc_v_T_B(v_s_g, v_s_l)
        v_L_L_S = self.calc_v_L_L_S(v_s_g, v_s_l)
        return v_T_B - (v_T_B - v_L_L_S)*H_L_L_S/H_L_T_B
    
    def calc_v_g_T_B(self,v_s_g, v_s_l):         #Actual Taylor Bubble velocity in Taylor Bubble Region
        v_s = self.calc_v_s(v_s_g, v_s_l)
        H_L_T_B = self.calc_H_L_T_B(self.h_f_d)
        v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        return (v_s - v_L_T_B*H_L_T_B)/(1 - H_L_T_B)
    
    #Геометрические параметры в для региона Тейлор бабл
    def calc_S_f(self,h_f_d):
        b = 2*h_f_d-1
        return self.d_pipe*(math.pi - math.acos(b))
    
    def calc_S_g(self,h_f_d):
        S_f = self.calc_S_f(h_f_d)
        return math.pi*self.d_pipe - S_f
    
    def calc_S_i(self,h_f_d):
        b = 2*h_f_d-1
        return self.d_pipe*(1-b*b)**0.5
    
    def calc_Area_film(self,h_f_d):
        Area_pipe = self.calc_Area_pipe(self.d_pipe)
        H_L_T_B = self.calc_H_L_T_B(h_f_d)
        return H_L_T_B*Area_pipe
    
    def calc_Area_gas(self,h_f_d):
        Area_pipe = self.calc_Area_pipe(self.d_pipe)
        H_L_T_B = self.calc_H_L_T_B(h_f_d)
        return (1-H_L_T_B)*Area_pipe   
    
    def calc_d_hydr_film(self,h_f_d):         #расчет гидравлического диаметра для пленки
        Area_film = self.calc_Area_film(h_f_d)
        S_f = self.calc_S_f(h_f_d)
        return 4*Area_film/S_f
    
    def calc_d_hydr_gas_bubble(self,h_f_d):         #расчет гидравлического диаметра для пузырька Тейлора
        Area_gas = self.calc_Area_gas(h_f_d)
        S_g = self.calc_S_g(h_f_d)
        S_i = self.calc_S_i(h_f_d)
        return 4*Area_gas/(S_g + S_i) 
    
    def calc_N_Re_film(self,v_s_g, v_s_l, h_f_d):     # Расчет числа Рейнольдса для пленки
        d_hydr_film = self.calc_d_hydr_film(h_f_d)
        v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        return self.rho_l_PT*d_hydr_film*math.fabs(v_L_T_B)/self.mu_l_PT
    
    def calc_N_Re_gas_bubble(self,v_s_g, v_s_l, h_f_d):     # Расчет числа Рейнольдса для пузырька Тейлора
        d_hydr_gas_bubble = self.calc_d_hydr_gas_bubble(h_f_d)
        v_g_T_B = self.calc_v_g_T_B(v_s_g, v_s_l)
        return self.rho_g_PT*d_hydr_gas_bubble*math.fabs(v_g_T_B)/self.mu_g_PT
     
    def calc_f_factor_film(self,v_s_g, v_s_l):     # Расчет гидравлического коэффициента потерь на трение для пленки
        if self.calc_N_Re_film(v_s_g, v_s_l, self.h_f_d) < 2000:
        # режим ламинарный
            return 64*self.calc_N_Re_film(v_s_g, v_s_l, self.h_f_d)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_film(v_s_g, v_s_l, self.h_f_d)**(-0.2)
     
    def calc_f_factor_gas_bubble(self,v_s_g, v_s_l):     # Расчет гидравлического коэффициента потерь на трение для пузырька Тейлора
        if self.calc_N_Re_gas_bubble(self.v_s_g, self.v_s_l, self.h_f_d)<2000:
        # режим ламинарный
            return 64*self.calc_N_Re_gas_bubble(self.v_s_g, self.v_s_l, self.h_f_d)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_gas_bubble(self.v_s_g, self.v_s_l, self.h_f_d)**(-0.2)
            
    #Shear Stress
    def calc_tau_w_film(self,v_s_g, v_s_l):      #Liquid Wall Shear Stress for film in Taylor Bubble Region
        f_factor_film = self.calc_f_factor_film(v_s_g, v_s_l) 
        v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        return f_factor_film*self.rho_l_PT*v_L_T_B*math.fabs(v_L_T_B)/8
    
    def calc_tau_w_gas_bubble(self,v_s_g, v_s_l):      #Gas Wall Shear Stress for gas bubble in Taylor Bubble Region
        f_factor_gas_bubble = self.calc_f_factor_gas_bubble(v_s_g, v_s_l)
        v_g_T_B = self.calc_v_g_T_B(v_s_g, v_s_l)
        return f_factor_gas_bubble*self.rho_g_PT*v_g_T_B*math.fabs(v_g_T_B)/8
    
    def calc_tau_i(self,v_s_g, v_s_l):       #Interfacial Shear Stress
        f_factor_i = 0.0568
        v_g_T_B = self.calc_v_g_T_B(v_s_g, v_s_l)
        v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        return f_factor_i*self.rho_g_PT*(v_g_T_B - v_L_T_B)*math.fabs(v_g_T_B - v_L_T_B)/8

    def combined_momentum_eq_film_reg(self,v_s_g, v_s_l):       #Combined momentum equation in film region
        self.h_f_d = 0.001
        res = 0
        while (res<0):
#        h_f_d=h_f_d+0.0001 
            tau_w_film= self.calc_tau_w_film(v_s_g, v_s_l)
            tau_w_gas_bubble = self.calc_tau_w_gas_bubble(v_s_g, v_s_l)
            tau_i = self.calc_tau_i(v_s_g, v_s_l)
            S_f = self.calc_S_f(self.h_f_d)
            S_g = self.calc_S_g(self.h_f_d)
            S_i = self.calc_S_i(self.h_f_d)
            Area_film = self.calc_Area_film(self.h_f_d)
            Area_gas = self.calc_Area_gas(self.h_f_d)
            a1 = tau_w_gas_bubble*S_g/Area_gas
            b1 = tau_w_film*S_f/Area_film
            c11 = tau_i*S_i
            c12 = 1/Area_gas + 1/Area_film
            c1 = c11*c12
            d1 = (self.rho_l_PT - self.rho_g_PT)*self.gravity*math.sin(math.pi*self.angle/180)
            #return b1
            res= a1-b1+c1-d1
            self.h_f_d= self.h_f_d+0.0001 
        return res
    
    def calc_H_L_T_B_fact(self,v_s_g, v_s_l):
        if self.h_f_d>=1:
            return 0.999
        elif self.h_f_d<=0:
            return 0.001
        else:
            return self.h_f_d      

    
    #h_f_d = fsolve(combined_momentum_eq_film_reg, 0.001)
    #h_f_d = 0.5
        
    
    
    #print(h_f_d)
        
    

    
    #print(calc_H_L_T_B_fact(h_f_d))
           
    #Определение Liquid Holdup в замисимости от структуры потока
    def calc_regime_liquid_holdup(self,v_s_g, v_s_l, rho_g_PT):
        
#        ff = self.calc_uniq_func_dim_height
        self.dim_h_l = fsolve(self.calc_uniq_func_dim_height, 0.001)

        
        flow_structure = self.calc_flow_structure(self.dim_h_l)
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        H_L_total_an = self.calc_H_L_gas_core(v_s_g, v_s_l)     #Без учета холдапа в пленке в аннулар режиме
    #    H_L_total_an = calc_H_L_total_an(v_s_g, v_s_l)
        H_L_dispers = self.calc_H_L_dispers(v_s_g, v_s_l)
        H_L_strat = self.calc_H_L_strat(self.dim_h_l)
        H_L_S_U = self.calc_H_L_S_U(v_s_g, v_s_l)
        H_L_T_B = self.calc_H_L_T_B_fact(self.v_s_g, self.v_s_l)
        d = [H_L_T_B,H_L_S_U,H_L_L_S]
        if flow_structure == 1:
            return H_L_total_an
        elif flow_structure == 2:
            return H_L_strat
        elif flow_structure == 3:
            return H_L_strat
        elif flow_structure == 4:
            return H_L_dispers
        else:
            return d
    
#    h_f_d = 0.01
#    while (self.combined_momentum_eq_film_reg(v_s_g, v_s_l, h_f_d)>0):
#        h_f_d=h_f_d+0.0001    
        
    #Модуль расчета длины и частоты пробки
    def calc_L_S_feet(self,d_pipe):
        d_pipe_inch = d_pipe/0.0254
        a = math.log(d_pipe_inch) 
        b = -25.4 + 28.5*a**0.1
        return math.exp(b)
    
    def calc_L_S_feet_max(self,d_pipe):
        L_S_feet = self.calc_L_S_feet(d_pipe)
        b = math.log(L_S_feet)
        c = 0.5*3.08 + b
        return math.exp(c)
    
    def calc_L_S(self,d_pipe): #Average Slug Length in meter
        L_S_feet = self.calc_L_S_feet(d_pipe)
        return 0.3048*L_S_feet
    
    def calc_L_S_max(self,d_pipe): #Maximum Slug Length in meter
        L_S_feet_max = self.calc_L_S_feet_max(d_pipe)
        return 0.3048*L_S_feet_max
    
    def calc_L_U(self,v_s_g, v_s_l, d_pipe):         #Slug Unit Length
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        H_L_T_B = self.calc_H_L_T_B_fact(v_s_g, v_s_l)
        v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        v_L_L_S = self.calc_v_L_L_S(v_s_g, v_s_l)
        L_S = self.calc_L_S(d_pipe)
        a = L_S*(v_L_L_S*H_L_L_S - v_L_T_B*H_L_T_B)
        b = v_s_l - v_L_T_B*H_L_T_B
        return a/b
    
    def calc_L_f(self,v_s_g, v_s_l, d_pipe):     #Film Length
        L_S = self.calc_L_S(d_pipe)
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        return L_U - L_S
    
    def calc_Freq_S(self,v_s_g, v_s_l, d_pipe):       #Slug Frequency
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        v_T_B = self.calc_v_T_B(v_s_g, v_s_l)
        return v_T_B/L_U
                       
    def calc_Freq_S_Zabaras(self,v_s_g, v_s_l, d_pipe):      #Empirical correlation for Slug Frequency by Zabaras (2000)     
        v_s_l_feet = v_s_l/0.3048
        gravity_feet = 32.2
        d_pipe_feet = d_pipe/0.3048
        v_mix_feet = self.calc_v_s_mix(v_s_l, v_s_g)/0.3048
        a = (v_s_l_feet/(gravity_feet*d_pipe_feet))**1.2
        b = (212.6/v_mix_feet + v_mix_feet)**1.2
        c = math.sin((math.pi*self.angle/180))
        d = 0.836 + 2.75*c**0.25
        return 0.0226*a*b*d
    
    #Модуль расчета градиента давления
    
    #Intermittent Flow (Slug Regime)
    def calc_pressure_drop_film_reg(self,v_s_g, v_s_l, d_pipe):      #Pressure Drop in Film Region
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        L_f = self.calc_L_f(v_s_g, v_s_l, d_pipe)
        tau_w_film= self.calc_tau_w_film(v_s_g, v_s_l)
        tau_w_gas_bubble = self.calc_tau_w_gas_bubble(v_s_g, v_s_l)
        S_f = self.calc_S_f(self.h_f_d)
        S_g = self.calc_S_g(self.h_f_d)
        Area_pipe = self.calc_Area_pipe(d_pipe)
        H_L_T_B = self.calc_H_L_T_B_fact(v_s_g, v_s_l)
        rho_f = self.rho_l_PT*H_L_T_B + self.rho_g_PT*(1-H_L_T_B) 
        a = L_f/L_U
        b = rho_f*self.gravity*math.sin((math.pi*self.angle/180))
        c = (tau_w_film*S_f + tau_w_gas_bubble*S_g)/Area_pipe
        d = a*(c + b)
        return -d
    
    def calc_mu_slug_PT(self,mu_l_PT, mu_g_PT, v_s_g, v_s_l):      #Вязкость слага
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        return mu_l_PT*H_L_L_S + mu_g_PT*(1-H_L_L_S)
    
    def calc_rho_slug_PT(self,v_s_g, v_s_l):         #Плотность слага
        H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        return self.rho_l_PT*H_L_L_S + self.rho_g_PT*(1-H_L_L_S) 
        
    def calc_N_Re_Slug(self,v_s_g, v_s_l):     # Расчет числа Рейнольдса для Slug
        v_s_mix = self.calc_v_s_mix(v_s_l, v_s_g)
        rho_slug = self.calc_rho_slug_PT(v_s_g, v_s_l)
        mu_slug_PT = self.calc_mu_slug_PT(self.mu_l_PT, self.mu_g_PT, v_s_g, v_s_l)
        return rho_slug*self.d_pipe*v_s_mix/mu_slug_PT
    
    def calc_f_factor_Slug(self,v_s_g, v_s_l):     # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.calc_N_Re_Slug(v_s_g, v_s_l)<2000:
        # режим ламинарный
            return 64*self.calc_N_Re_Slug(v_s_g, v_s_l)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_N_Re_Slug(v_s_g, v_s_l)**(-0.2)
        
    def calc_tau_Slug(self,v_s_g, v_s_l):        #Shear Stress in Slug Region
        f_factor_Slug = self.calc_f_factor_Slug(v_s_g, v_s_l)
        rho_slug = self.calc_rho_slug_PT(v_s_g, v_s_l)
        v_s_mix = self.calc_v_s_mix(v_s_l, v_s_g)
        return f_factor_Slug*rho_slug*v_s_mix*v_s_mix/8
    
    def calc_pressure_drop_slug_reg(self,v_s_g, v_s_l, d_pipe):      #Pressure Drop in Slug Region
        tau_Slug = self.calc_tau_Slug(v_s_g, v_s_l)
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        L_S = self.calc_L_S(d_pipe)
        rho_slug = self.calc_rho_slug_PT(v_s_g, v_s_l)
        a = L_S/L_U
        b = tau_Slug*4/self.d_pipe
        c = rho_slug*self.gravity*math.sin((math.pi*self.angle/180))
        d = a*(b + c)
        return -d
    
    def calc_pressure_drop_S_U(self,v_s_g, v_s_l, d_pipe):        #Pressure Drop in Intermittent Flow
        a = self.calc_pressure_drop_slug_reg(v_s_g, v_s_l, d_pipe)
        b = self.calc_pressure_drop_film_reg(v_s_g, v_s_l, d_pipe)
        return a+b
    
    #Dispersed Bubble Flow
    def calc_rho_no_slip(self,v_s_g, v_s_l):
        H_L_dispers = self.calc_H_L_dispers(v_s_g, v_s_l)
        return self.rho_l_PT*H_L_dispers + self.rho_g_PT*(1-H_L_dispers)
    
    def calc_mu_no_slip(self,mu_l_PT, mu_g_PT):
        H_L_dispers = self.calc_H_L_dispers(self.v_s_g, self.v_s_l) 
        return H_L_dispers*mu_l_PT + (1-H_L_dispers)*mu_g_PT
    
    def calc_Re_no_slip(self,v_s_g, v_s_l, d_pipe):
        rho_no_slip = self.calc_rho_no_slip(v_s_g, v_s_l)
        mu_no_slip = self.calc_mu_no_slip(self.mu_l_PT, self.mu_g_PT)
        v_s_mix = self.calc_v_s_mix(v_s_l, v_s_g)
        return v_s_mix*d_pipe*rho_no_slip/mu_no_slip
    
    def calc_f_factor_no_slip(self,v_s_g, v_s_l):     # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.calc_Re_no_slip(v_s_g, v_s_l, self.d_pipe)<2000:
        # режим ламинарный
            return 64*self.calc_Re_no_slip(v_s_g, v_s_l, self.d_pipe)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_Re_no_slip(v_s_g, v_s_l, self.d_pipe)**(-0.2)
    
    def calc_tau_no_slip(self,v_s_g, v_s_l):        #Shear Stress in Dispersed Region
        f_factor_no_slip = self.calc_f_factor_no_slip(v_s_g, v_s_l)
        rho_no_slip = self.calc_rho_no_slip(v_s_g, v_s_l)
        v_s_mix = self.calc_v_s_mix(v_s_l, v_s_g)
        return f_factor_no_slip*rho_no_slip*v_s_mix*v_s_mix/8
    
    def calc_pressure_drop_dispers(self,v_s_g, v_s_l, d_pipe):        #Pressure Drop in Dispersed Bubble Flow
        rho_no_slip = self.calc_rho_no_slip(v_s_g, v_s_l)
        tau_no_slip = self.calc_tau_no_slip(v_s_g, v_s_l)
        a = tau_no_slip*4/d_pipe
        b = rho_no_slip*self.gravity*math.sin((math.pi*self.angle/180))
        c = a+b
        return -c
    
    #Stratified Flow (Smooth and Wavy)
    def calc_S_L_Str(self,dim_h_l):
        h_l=self.calc_dim_h_l_fact(dim_h_l)
        b = 2*h_l-1
        return self.d_pipe*(math.pi - math.acos(b))
    
    def calc_S_g_Str(self,dim_h_l):
        S_L_Str = self.calc_S_L_Str(dim_h_l)
        return math.pi*self.d_pipe - S_L_Str
    
    def calc_S_i_Str(self,dim_h_l):
        h_l=self.calc_dim_h_l_fact(dim_h_l)
        b = 2*h_l-1
        return self.d_pipe*(1-b*b)**0.5
    
    def calc_Area_liquid_Str(self,dim_h_l):
        Area_pipe = self.calc_Area_pipe(self.d_pipe)
        H_L_strat = self.calc_H_L_strat(dim_h_l)
        return H_L_strat*Area_pipe
    
    def calc_Area_gas_Str(self,dim_h_l):
        Area_pipe = self.calc_Area_pipe(self.d_pipe)
        H_L_strat = self.calc_H_L_strat(dim_h_l)
        return (1-H_L_strat)*Area_pipe   
    
    def calc_d_hydr_liq(self,dim_h_l):         #расчет гидравлического диаметра для слоя жидкости
        Area_liquid_Str = self.calc_Area_liquid_Str(dim_h_l)
        S_L_Str = self.calc_S_L_Str(dim_h_l)
        return 4*Area_liquid_Str/S_L_Str
    
    def calc_d_hydr_gas(self,dim_h_l):         #расчет гидравлического диаметра для слоя газа
        Area_gas = self.calc_Area_gas_Str(dim_h_l)
        S_g = self.calc_S_g_Str(dim_h_l)
        S_i = self.calc_S_i_Str(dim_h_l)
        return 4*Area_gas/(S_g + S_i) 
    
    def calc_v_L_Str(self,v_s_l):         #Actual liquid velocity
        H_L_strat = self.calc_H_L_strat(self.dim_h_l)
        return v_s_l/H_L_strat
    
    def calc_v_g_Str(self,v_s_g):         #Actual gas velocity
        H_L_strat = self.calc_H_L_strat(self.dim_h_l)
        return v_s_g/(1-H_L_strat)
    
    def calc_Re_gas_Str(self,v_s_g): 
        v_g_Str = self.calc_v_g_Str(v_s_g)
        d_hydr_gas = self.calc_d_hydr_gas(self.dim_h_l)
        return v_g_Str*d_hydr_gas*self.rho_g_PT/self.mu_g_PT
    
    def calc_Re_liquid_Str(self,v_s_l): 
        v_L_Str = self.calc_v_L_Str(v_s_l)
        d_hydr_liq = self.calc_d_hydr_liq(self.dim_h_l)
        return v_L_Str*d_hydr_liq*self.rho_l_PT/self.mu_l_PT
    
    def calc_f_factor_l_str(self,v_s_l):     # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.calc_Re_liquid_Str(v_s_l)<2000:
        # режим ламинарный
            return 64*self.calc_Re_liquid_Str(v_s_l)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_Re_liquid_Str(v_s_l)**(-0.2)
    
    def calc_f_factor_g_str(self,v_s_g):     # Расчет гидравлического коэффициента потерь на трение для газа
        if self.calc_Re_gas_Str(v_s_g)<2000:
        # режим ламинарный
            return 64*self.calc_Re_gas_Str(v_s_g)**(-1)
        else:
            # режим турбулентный
            return 0.184*self.calc_Re_gas_Str(v_s_g)**(-0.2)
    
    def calc_f_factor_i_str(self,dim_h_l):
        a = self.calc_flow_structure(dim_h_l)
        b = self.calc_f_factor_g_str(self.v_s_g)
        c = 0.0568
        if a == 2:
            return c
        elif a == 3:
            return b
        else:
            return 0.2
    
    
    def calc_tau_w_L_str(self,v_s_l):
        f_factor_l_str = self.calc_f_factor_l_str(v_s_l)
        v_L_Str = self.calc_v_L_Str(v_s_l)
        return f_factor_l_str*self.rho_l_PT*v_L_Str*v_L_Str/8
    
    def calc_tau_w_g_str(self,v_s_g):
        f_factor_g_str = self.calc_f_factor_g_str(v_s_g)
        v_g_Str = self.calc_v_g_Str(v_s_g)
        return f_factor_g_str*self.rho_g_PT*v_g_Str*v_g_Str/8
    
    def calc_tau_i_str(self,v_s_l, v_s_g):
        f_factor_i_str = self.calc_f_factor_i_str(self.dim_h_l)
        v_L_Str = self.calc_v_L_Str(v_s_l)
        v_g_Str = self.calc_v_g_Str(v_s_g)
        return 0.125*f_factor_i_str*self.rho_g_PT*(v_g_Str-v_L_Str)**2
    
    def calc_pressure_drop_gas_str(self,v_s_g, v_s_l, d_pipe):
        tau_i_str = self.calc_tau_i_str(v_s_l, v_s_g)
        tau_w_g_str = self.calc_tau_w_g_str(v_s_g)
        S_i = self.calc_S_i_Str(self.dim_h_l)
        S_g = self.calc_S_g_Str(self.dim_h_l)
        Area_gas = self.calc_Area_gas_Str(self.dim_h_l)
        a = tau_i_str*S_i/Area_gas
        b = tau_w_g_str*S_g/Area_gas
        c = self.rho_g_PT*self.gravity*math.sin((math.pi*self.angle/180))
        d = a + b + c
        return -d
    
    def calc_pressure_drop_liquid_str(self,v_s_g, v_s_l, d_pipe):
        tau_i_str = self.calc_tau_i_str(v_s_l, v_s_g)
        tau_w_L_str = self.calc_tau_w_L_str(v_s_l)
        S_i = self.calc_S_i_Str(self.dim_h_l)
        S_L = self.calc_S_L_Str(self.dim_h_l)
        Area_liquid_Str = self.calc_Area_liquid_Str(self.dim_h_l)
        a = -tau_i_str*S_i/Area_liquid_Str
        b = tau_w_L_str*S_L/Area_liquid_Str
        c = self.rho_l_PT*self.gravity*math.sin((math.pi*self.angle/180))
        d = a + b + c
        return -d
    
    def calc_pressure_drop_str(self,v_s_g, v_s_l, d_pipe):
        a = self.calc_pressure_drop_liquid_str(v_s_g, v_s_l, d_pipe)
        b = self.calc_pressure_drop_gas_str(v_s_g, v_s_l, d_pipe)
        return a+b
    
    #Расчет градиета давления
    def pressure_drop(self,v_s_g, v_s_l, d_pipe):
        pressure_drop_str = self.calc_pressure_drop_str(v_s_g, v_s_l, d_pipe)
        pressure_drop_dispers = self.calc_pressure_drop_dispers(v_s_g, v_s_l, d_pipe)
        pressure_drop_S_U = self.calc_pressure_drop_S_U(v_s_g, v_s_l, d_pipe)
        a = self.calc_flow_structure(self.dim_h_l)
        if a == 2:
            return pressure_drop_str
        elif a == 3:
            return pressure_drop_str
        elif a == 4:
            return pressure_drop_dispers
        elif a == 5:
            return pressure_drop_S_U
    """    
            # приведенная скорость жидкости, м3/с
                # приведенная скорость газа, м3/с
                # давление в трубе, МПа
             # температура на устье скважины, К
              # длина трубы, м
        # диаметр трубы, м
          # плотность жидкости при термобарических условиях, кг/м3
         # плотность газа при термобарических условиях, кг/м3
        # вязкость жидкости при термобарических условиях, Па*с
        # вязкость газа при термобарических условиях, Па*с
         # угол наклона трубы, градус
       # абсолютная шероховатость трубы, м
     #ускорение свободного падения, м/с2
    #поверхностное натяжение на границе жидкости и газа, Н/м    
    """    

grX = grad_Xiao()

dP = 0.4
print (grX.calc_L_S(dP))
print (grX.calc_L_S_max(dP))
a = grX.combined_momentum_eq_film_reg(1,1)
b= grX.calc_regime_liquid_holdup(1,1,0.6)
print (grX.calc_L_U(2, 10, dP))
print (grX.pressure_drop(1, 1, 0.4))

"""
print (grX.calc_Freq_S(v_s_g, v_s_l, d_pipe))
print (grX.calc_Freq_S_Zabaras(v_s_g, v_s_l, d_pipe))
print (grX.calc_flow_structure(dim_h_l))
print (grX.calc_pressure_drop_S_U(v_s_g, v_s_l, d_pipe))
print (grX.pressure_drop(1, 1, 0.4))

t=np.arange(0.15, 0.95, 0.001)
yy=[]

dim_h_l = fsolve(calc_uniq_func_dim_height, 0.001)


vsl=0.5
d=0.3

#print (h_f_d)
t=np.arange(0.1, 10, 0.01)
yy=[]
 
for vsg in t :
    yy.append(self.pressure_drop(vsg, vsl, d, angle=90))
plt.plot(t,yy) 
plt.grid(True)


"""