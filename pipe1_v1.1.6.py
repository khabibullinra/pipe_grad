# -*- coding: utf-8 -*-
"""
Редактор Spyder
 
Это временный скриптовый файл.

пример изменений

"""

# Загрузка библиотек необходимых для отрисовки графиков


import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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


class GradXiao:

    # Расчет предварительных величин
    def __init__(self):
        """
        конструктор класса, задает исходные значения 
        """
        # Исходные данные для проведения расчета
        self.v_s_l = 0.2          # приведенная скорость жидкости, м3/с
        self.v_s_g = 10           # приведенная скорость газа, м3/с
        self.P_pipe = 40            # давление в трубе, МПа
        self.T_pipe = 350          # температура на устье скважины, К
        self.L_pipe = 2000          # длина трубы, м
        self.d_pipe = 0.8       # диаметр трубы, м
        self.rho_l_PT = 1000      # плотность жидкости при термобарических условиях, кг/м3
        self.rho_g_PT = 20      # плотность газа при термобарических условиях, кг/м3
        self.mu_l_PT = 0.001    # вязкость жидкости при термобарических условиях, Па*с
        self.mu_g_PT = 0.00002    # вязкость газа при термобарических условиях, Па*с
        self.angle = 0          # угол наклона трубы, градус
        self.rough = 0.0001     # абсолютная шероховатость трубы, м
        self.gravity = 9.81     # ускорение свободного падения, м/с2
        self.sigma_l = 0.03     # поверхностное натяжение на границе жидкости и газа, Н/м

    def calc_f_factor_l(self):     # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.N_Re_l < 2000:
            """
            режим ламинарный
            """
            return 64*self.N_Re_l**(-1)
        else:
            """
            режим турбулентный
            """
            return 0.184*self.N_Re_l**(-0.2)
     
    def calc_f_factor_g(self):     # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.N_Re_g<2000:
            # режим ламинарный
            return 64*self.N_Re_g**(-1)
        else:
            # режим турбулентный
            return 0.184*self.N_Re_g**(-0.2)
        
    def calc_auxiliary_init_vars(self):  # расчет вспомогательных переменных
        self.Area_pipe_m =  math.pi*self.d_pipe**2/4  # диаметр трубы
        self.d_pipe_inch = self.d_pipe/0.0254
        # расчет длины слага, зависит только от диаметра
        a = math.log(self.d_pipe_inch) 
        b = -25.4 + 28.5*a**0.1
        self.L_S_feet = math.exp(b)
        self.L_S = self.L_S_feet * 0.3048         
        # расчет максимальной длины слага 
        b = math.log(self.L_S_feet)
        c = 0.5*3.08 + b
        self.L_S_feet_max = math.exp(c) 
        self.L_S_max = self.L_S_feet_max * 0.3048         
        self.v_s_mix = self.v_s_g + self.v_s_l
        # найдем число Рейнольдса для жидкости для заданных условий
        self.N_Re_l = self.rho_l_PT * self.d_pipe*self.v_s_l / self.mu_l_PT
        self.N_Re_g = self.rho_g_PT * self.d_pipe * self.v_s_g / self.mu_g_PT
        self.N_Re_mix_slug = self.v_s_mix * self.rho_l_PT * self.d_pipe / self.mu_l_PT
        self.N_Re_mix_slug_crit = 1000
        self.c_o = 2.27/(1+(self.N_Re_mix_slug / self.N_Re_mix_slug_crit)**2) + 1.2 / (1+(self.N_Re_mix_slug_crit / self.N_Re_mix_slug)**2)
        aa = (self.gravity*self.d_pipe)**0.5
        self.v_d_TB = 0.54 * aa * math.cos(math.pi * self.angle/180) + 0.35 * aa ** 0.5 * math.sin(math.pi * self.angle / 180)
        self.v_T_B = self.c_o * self.v_s_mix + self.v_d_TB
        self.H_L_L_S = 1/(1+(self.v_s_mix/8.66)**1.39)
        self.H_L_dispers = self.v_s_l / (self.v_s_l + self.v_s_g)
        c = (self.gravity*self.sigma_l*(self.rho_l_PT-self.rho_g_PT)/self.rho_l_PT**2)**0.25
        self.v_g_L_S = self.v_s_mix + 1.53 * c * self.H_L_L_S**0.1 * math.sin(math.pi * self.angle / 180)
        self.v_L_L_S = (self.v_s_l + self.v_s_g - self.v_g_L_S * (1 - self.H_L_L_S)) / self.H_L_L_S
        self.H_L_S_U = (self.v_T_B*self.H_L_L_S + self.v_g_L_S * (1-self.H_L_L_S) - self.v_s_g) / self.v_T_B

        if self.N_Re_l<2000:
            self.n_power_coeff = 1
        else:
            self.n_power_coeff = 0.2       

        if self.N_Re_g<2000:
            self.m_power_coeff = 1
        else:
            self.m_power_coeff = 0.2    

        self.f_factor_l = self.calc_f_factor_l()
        self.f_factor_g = self.calc_f_factor_g()
        self.pres_drop_s_l = - self.f_factor_l * self.rho_l_PT * self.v_s_l**2 / (2 * self.d_pipe)
        self.pres_drop_s_g = - self.f_factor_g * self.rho_g_PT * self.v_s_g**2 / (2 * self.d_pipe)
        self.x_large = ((-self.pres_drop_s_l)/(-self.pres_drop_s_g))**0.5
        self.y_large = (self.rho_l_PT-self.rho_g_PT)*self.gravity*math.sin(math.pi*self.angle/180)/(-self.pres_drop_s_g)        
        self.f_large = (self.rho_g_PT / (self.rho_l_PT - self.rho_g_PT))**0.5 * self.v_s_g/(self.d_pipe * self.gravity * math.cos(math.pi * self.angle/180))**0.5
        self.k_large = self.N_Re_l**0.5*self.f_large
        self.t_large = (-self.pres_drop_s_l/((self.rho_l_PT-self.rho_g_PT)*self.gravity*math.cos(math.pi*self.angle/180)))**0.5
        #Slug Liquid Holdup in Slug Region (fig 4.14)
        return

    def calc_auxiliary_hfd_vars(self):
        """
        Здесь расчитываются параметры, которые надо пересчитать после оценнки толщины пленки
        :return:
        """

        # Определение Liquid Holdup in Taylor Bubble
        # Вычисление холдапа жидкоти в Тейлор бабл регионе из предположения размера пленки, как для стратифайд режима
        self.H_L_T_B = 0.5 + (4 / math.pi) * (self.h_f_d - 0.5) * (self.h_f_d - self.h_f_d ** 2) ** 0.5 + (
                                                                         1 / math.pi) * math.asin(2 * self.h_f_d - 1)
        # Actual film velocity in Taylor Bubble Region
        self.v_L_T_B = self.v_T_B - (self.v_T_B - self.v_L_L_S) * self.H_L_L_S / self.H_L_T_B
        # Actual Taylor Bubble velocity in Taylor Bubble Region
        self.v_g_T_B = (self.v_s_mix - self.v_L_T_B * self.H_L_T_B) / (1 - self.H_L_T_B)
        # Геометрические параметры в для региона Тейлор бабл
        b = 2 * self.h_f_d - 1
        self.S_f = self.d_pipe * (math.pi - math.acos(b))
        self.S_g = math.pi * self.d_pipe - self.S_f
        b = 2 * self.h_f_d - 1
        self.S_i = self.d_pipe * (1 - b * b) ** 0.5
        self.Area_film = self.H_L_T_B * self.Area_pipe_m
        self.Area_gas = (1 - self.H_L_T_B) * self.Area_pipe_m
        self.d_hydr_film = 4 * self.Area_film / self.S_f # расчет гидравлического диаметра для пленки
        self.d_hydr_gas_bubble = 4 * self.Area_gas / (self.S_g + self.S_i)# расчет гидравлического диаметра для пузырька Тейлора
        self.N_Re_film = self.rho_l_PT * self.d_hydr_film * math.fabs(self.v_L_T_B) / self.mu_l_PT# Расчет числа Рейнольдса для пленки
        self.N_Re_gas_bubble = self.rho_g_PT * self.d_hydr_gas_bubble * math.fabs(self.v_g_T_B) / self.mu_g_PT# Расчет числа Рейнольдса для пузырька Тейлора
        if self.N_Re_film < 2000:# Расчет гидравлического коэффициента потерь на трение для пленки
            # режим ламинарный
            self.f_factor_film = 64 * self.N_Re_film ** (-1)
        else:
            # режим турбулентный
            self.f_factor_film =  0.184 * self.N_Re_film ** (-0.2)
        #Расчет гидравлического коэффициента потерь на трение для пузырька Тейлора
        if self.N_Re_gas_bubble < 2000:
            # режим ламинарный
            self.f_factor_gas_bubble = 64 * self.N_Re_gas_bubble ** (-1)
        else:
            # режим турбулентный
            self.f_factor_gas_bubble =  0.184 * self.N_Re_gas_bubble ** (-0.2)
        # Shear Stress
        # Liquid Wall Shear Stress for film in Taylor Bubble Region
        self.tau_w_film = self.f_factor_film * self.rho_l_PT * self.v_L_T_B * math.fabs(self.v_L_T_B) / 8
        # def calc_tau_w_gas_bubble(self):
        #  Gas Wall Shear Stress for gas bubble in Taylor Bubble Region
        self.tau_w_gas_bubble = self.f_factor_gas_bubble * self.rho_g_PT * self.v_g_T_B * math.fabs(self.v_g_T_B) / 8
        #def calc_tau_i(self, v_s_g, v_s_l):
        # Interfacial Shear Stress
        f_factor_i = 0.0568
        self.tau_i = f_factor_i * self.rho_g_PT * (self.v_g_T_B - self.v_L_T_B) * math.fabs(self.v_g_T_B - self.v_L_T_B) / 8

        #def calc_L_U(self, v_s_g, v_s_l, d_pipe):  # Slug Unit Length
        #    H_L_L_S = self.calc_H_L_L_S(v_s_g, v_s_l)
        #    H_L_T_B = self.calc_H_L_T_B_fact(v_s_g, v_s_l)
        #    v_L_T_B = self.calc_v_L_T_B(v_s_g, v_s_l)
        #    v_L_L_S = self.calc_v_L_L_S(v_s_g, v_s_l)
        #    L_S = self.L_S
        a = self.L_S * (self.v_L_L_S * self.H_L_L_S - self.v_L_T_B * self.H_L_T_B)
        b = self.v_s_l - self.v_L_T_B * self.H_L_T_B
        self.L_U = a / b

        #def calc_L_f(self, v_s_g, v_s_l, d_pipe):  # Film Length
        #    L_S = self.L_S
        #    L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        self.L_f = self.L_U - self.L_S

        #def calc_Freq_S(self, v_s_g, v_s_l, d_pipe):  # Slug Frequency
        #    L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        #    # v_T_B = self.calc_v_T_B(v_s_g, v_s_l)
        self.Freq_S = self.v_T_B / self.L_U

        #def calc_Freq_S_Zabaras(self, v_s_g, v_s_l,
        #                        d_pipe):  # Empirical correlation for Slug Frequency by Zabaras (2000)
        v_s_l_feet = self.v_s_l / 0.3048
        gravity_feet = 32.2
        d_pipe_feet = self.d_pipe / 0.3048
        v_mix_feet = self.v_s_mix / 0.3048
        a = (v_s_l_feet / (gravity_feet * d_pipe_feet)) ** 1.2
        b = (212.6 / v_mix_feet + v_mix_feet) ** 1.2
        c = math.sin((math.pi * self.angle / 180))
        d = 0.836 + 2.75 * c ** 0.25
        self.S_Zabaras = 0.0226 * a * b * d

        return

    def calc_dim_h_l_fact(self, dim_h_l):
        if dim_h_l >= 1:
            return 0.9999
        elif dim_h_l <= 0:
            return 0.0001
        else:
            return dim_h_l

    def calc_auxiliary_dimhl_vars(self):
        # определение струтуры потока
        """
        Stratified to Nonstratified Transition criterea (если критерий больше либо равен единице, то Nonstratified)
        """
        trans_A = self.f_large ** 2 * (1 / (1 - self.dim_h_l) ** 2 * (self.dim_v_g ** 2 * self.dim_S_i / self.dim_A_g))
        # Intermittent or Dispersed_bubble to Annular Transition criterea (если критерий меньше либо равен 0,35, то Annular)
        trans_B = self.dim_h_l
        # Stratified-Smooth to Stratified-Wavy Transition criterea (если больше либо равно нулю, то Stratified-Smooth)
        s=0.01
        trans_C = self.k_large ** 2 - 2 / (self.dim_v_l ** 0.5 * self.dim_v_g * s ** 0.5)
        # Intermittent to Dispersed-Bubble Transition criterea (если критерий больше либо равен нулю, то Dispersed-Bubble)
        trans_D = self.t_large ** 2 - (8 * self.dim_A_g / (
            self.dim_S_i * self.dim_v_l ** 2 * (self.dim_v_l * self.dim_d_l) ** (-self.n_power_coeff)))
        """
        # если режим определен, что stratified, то выбираем wavy or smooth
        # 1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        """
        if trans_C < 0:
            stratified_regime = 3
        else:
            stratified_regime = 2
        # check dispersed-bubble or intermittent(i.e. slug or plug)
        #    # 1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        if trans_D < 0:
            bubble_slug_regime = 5
        else:
            bubble_slug_regime = 4
        # если режим определен, что nonstratified, то выбираем между annular or intermittent or dispersed-bubble
        # 1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
        if trans_A >= 1 and trans_B <= 0.35:
            #        return 1        #Annular Regime будет наблюдаться в трубах малого диаметра, где есть эффект капилляра
            nonstratified_regime = 2  # Мы принимаем, что у нас Stratified-Wavy Regime, в связи с тем что площадь стенки
            # очень большая жидкость из-за гравитации будет стекать, не наблюдаетяс эффект капилляра в большой трубе
        else:
            nonstratified_regime = bubble_slug_regime
            # 1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
            #    trans_A = self.calc_trans_A(dim_h_l)
        if trans_A < 1:
            self.flow_structure = stratified_regime
        else:
            self.flow_structure = nonstratified_regime

        # Stratified flow
        self.H_L_strat = 0.5 + (4 / math.pi) * (self.dim_h_l - 0.5) * (self.dim_h_l - self.dim_h_l ** 2) ** 0.5 + (1
                                                                   / math.pi) * math.asin(2 * self.dim_h_l - 1)

        return  self.flow_structure

    def calc_auxiliary_pres_drop(self):
        #Slug Flow
        H_L_T_B = self.calc_H_L_T_B_fact()
        self.rho_f = self.rho_l_PT*H_L_T_B + self.rho_g_PT*(1-H_L_T_B)
        self.mu_slug_PT = self.mu_l_PT*self.H_L_L_S + self.mu_g_PT*(1-self.H_L_L_S)      #Вязкость слага
        self.rho_slug_PT =  self.rho_l_PT*self.H_L_L_S + self.rho_g_PT*(1-self.H_L_L_S)        #Плотность слага
        self.N_Re_Slug = self.rho_slug_PT*self.d_pipe*self.v_s_mix/self.mu_slug_PT     # Расчет числа Рейнольдса для Slug
        # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.N_Re_Slug<2000:
        # режим ламинарный
            self.f_factor_Slug = 64*self.N_Re_Slug**(-1)
        else:
            # режим турбулентный
            self.f_factor_Slug = 0.184*self.N_Re_Slug**(-0.2)
            
        self.tau_Slug = self.f_factor_Slug*self.rho_slug_PT*self.v_s_mix*self.v_s_mix/8 #Shear Stress in Slug Region
          
        #Dispersed Bubble Flow
        self.rho_no_slip = self.rho_l_PT * self.H_L_dispers + self.rho_g_PT*(1 - self.H_L_dispers)
        self.mu_no_slip = self.H_L_dispers*self.mu_l_PT + (1-self.H_L_dispers)*self.mu_g_PT
        self.Re_no_slip = self.v_s_mix*self.d_pipe*self.rho_no_slip/self.mu_no_slip
        # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.Re_no_slip<2000:
        # режим ламинарный
            self.f_factor_no_slip = 64*self.Re_no_slip**(-1)
        else:
            # режим турбулентный
            self.f_factor_no_slip = 0.184*self.Re_no_slip**(-0.2)
        
        self.tau_no_slip = self.f_factor_no_slip*self.rho_no_slip*self.v_s_mix*self.v_s_mix/8        #Shear Stress in Dispersed Region
        
        #Stratified Flow (Smooth and Wavy)
        h_l=self.calc_dim_h_l_fact(self.dim_h_l)
        b = 2*h_l-1
        self.S_L_Str = self.d_pipe*(math.pi - math.acos(b))
        self.S_g_Str = math.pi*self.d_pipe - self.S_L_Str
        self.S_i_Str = self.d_pipe*(1-b*b)**0.5
        self.Area_liquid_Str = self.H_L_strat*self.Area_pipe_m
        self.Area_gas_Str = (1-self.H_L_strat)*self.Area_pipe_m
        self.d_hydr_liq = 4*self.Area_liquid_Str/self.S_L_Str         #ПРОВЕРИТЬ! расчет гидравлического диаметра для слоя жидкости
        self.d_hydr_gas = 4*self.Area_gas/(self.S_g_Str + self.S_i_Str) #ПРОВЕРИТЬ #расчет гидравлического диаметра для слоя газа
    
        self.v_L_Str = self.v_s_l/self.H_L_strat       #Actual liquid velocity
        self.v_g_Str = self.v_s_g/(1-self.H_L_strat)   ##Actual gas velocity
    
        self.N_Re_gas_Str = self.v_g_Str*self.d_hydr_gas*self.rho_g_PT/self.mu_g_PT
        self.N_Re_liquid_Str = self.v_L_Str*self.d_hydr_liq*self.rho_l_PT/self.mu_l_PT

             # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.N_Re_liquid_Str<2000:
        # режим ламинарный
            self.f_factor_l_str = 64*self.N_Re_liquid_Str**(-1)
        else:
            # режим турбулентный
            self.f_factor_l_str = 0.184*self.N_Re_liquid_Str**(-0.2)
    
            # Расчет гидравлического коэффициента потерь на трение для газа
        if self.N_Re_gas_Str<2000:
        # режим ламинарный
            self.f_factor_g_str = 64*self.N_Re_gas_Str**(-1)
        else:
            # режим турбулентный
            self.f_factor_g_str = 0.184*self.N_Re_gas_Str**(-0.2)
 
    # Расчет гидравлического коэффициента потерь на трение между жидкостью и газом
        aa = self.flow_structure
        cc = 0.0568
        if aa == 2:
            self.f_factor_i_str = cc
        elif aa == 3:
            self.f_factor_i_str = self.f_factor_g_str
        else:
            self.f_factor_i_str = 0

        self.tau_w_L_str = self.f_factor_l_str*self.rho_l_PT*self.v_L_Str*self.v_L_Str/8
        self.tau_w_g_str = self.f_factor_g_str*self.rho_g_PT*self.v_g_Str*self.v_g_Str/8
        self.tau_i_str = 0.125*self.f_factor_i_str*self.rho_g_PT*(self.v_g_Str-self.v_L_Str)*(self.v_g_Str-self.v_L_Str)
       
        
        return

    def calc_uniq_func_dim_height(self,dim_h_l):
        def calc_dim_A_l(dim_h_l):
            b=2*dim_h_l-1
            a=math.pi - math.acos(b)
            c=1-b*b
            self.dim_A_l =  0.25*(a+b*c**0.5)
            return self.dim_A_l
         
        def calc_dim_A_g(dim_h_l):
            b=2*dim_h_l-1
            a=math.acos(b)
            c=1-b**2
            self.dim_A_g = 0.25*(a-b*c**0.5)
            return self.dim_A_g
         
        def calc_dim_S_l(dim_h_l):
            self.dim_S_l = math.pi - math.acos(2*dim_h_l-1)
            return self.dim_S_l
         
        def calc_dim_S_i(dim_h_l):
            self.dim_S_i = (1-(2*dim_h_l-1)**2)**0.5
            return self.dim_S_i 
           
        def calc_dim_S_g(dim_h_l):
            self.dim_S_g = math.acos(2*dim_h_l-1)
            return self.dim_S_g
           
        def calc_dim_v_l(dim_h_l):
            dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            self.dim_v_l = (dim_A_g/dim_A_l+1)
            return self.dim_v_l
           
        def calc_dim_v_g(dim_h_l):
            dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            self.dim_v_g = 1+dim_A_l/dim_A_g
            return self.dim_v_g 

        def calc_dim_d_l(dim_h_l):
            dim_A_l=0.25*(math.pi - math.acos(2*dim_h_l-1)+(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            dim_S_l=math.pi - math.acos(2*dim_h_l-1)
            self.dim_d_l = 4*dim_A_l/dim_S_l
            return self.dim_d_l
           
        def calc_dim_d_g(dim_h_l):
            dim_A_g=0.25*(math.acos(2*dim_h_l-1)-(2*dim_h_l-1)*(1-(2*dim_h_l-1)**2)**0.5)
            dim_S_g=math.acos(2*dim_h_l-1)
            dim_S_i=(1-(2*dim_h_l-1)**2)**0.5
            self.dim_d_g =  4*dim_A_g/(dim_S_i+dim_S_g) 
            return self.dim_d_g
        
        a1=self.x_large
        a2=self.y_large
        
        dim_h_l=self.calc_dim_h_l_fact(dim_h_l)
        b1=calc_dim_A_l(dim_h_l)
        b2=calc_dim_A_g(dim_h_l)
        c1=calc_dim_S_l(dim_h_l)
        c2=calc_dim_S_i(dim_h_l)
        c3=calc_dim_S_g(dim_h_l)
        d1=calc_dim_v_l(dim_h_l)
        d2=calc_dim_v_g(dim_h_l)
        e1=calc_dim_d_l(dim_h_l)
        e2=calc_dim_d_g(dim_h_l)
        f1=self.n_power_coeff
        f2=self.m_power_coeff
        a11=a1*a1
        b11=d1*e1
        c11=b11**(-f1)
        d11=d1*d1
        e11=c1/b1
        f11=d2*e2
        g11=f11**(-f2)
        h11=d2*d2
        i11=c3/b2+c2/b1+c2/b2
        res = a11*c11*d11*e11-g11*h11*i11+4*a2
        self.dim_h_l=dim_h_l
        return  res



    def combined_momentum_eq_film_reg(self):       #Combined momentum equation in film region
        self.h_f_d = 0.001
        res = 0.1
        while (res>0):
#        h_f_d=h_f_d+0.0001
            self.calc_auxiliary_hfd_vars()
            tau_w_film= self.tau_w_film
            tau_w_gas_bubble = self.tau_w_gas_bubble
            tau_i = self.tau_i
            S_f = self.S_f
            S_g = self.S_g
            S_i = self.S_i
            Area_film = self.Area_film
            Area_gas = self.Area_gas
            a1 = tau_w_gas_bubble*S_g/Area_gas
            b1 = tau_w_film*S_f/Area_film
            c11 = tau_i*S_i
            c12 = 1/Area_gas + 1/Area_film
            c1 = c11*c12
            d1 = (self.rho_l_PT - self.rho_g_PT)*self.gravity*math.sin(math.pi*self.angle/180)
            #return b1
            res= a1-b1+c1-d1
            self.h_f_d= self.h_f_d+0.00001 
        return res
    
    def calc_H_L_T_B_fact(self):
        if self.h_f_d>=1:
            return 0.999
        elif self.h_f_d<=0:
            return 0.001
        else:
            return self.h_f_d
           
    #Определение Liquid Holdup в замисимости от структуры потока
    def calc_regime_liquid_holdup(self):
        
        self.dim_h_l = fsolve(self.calc_uniq_func_dim_height, 0.001)

        self.calc_auxiliary_dimhl_vars()

        self.combined_momentum_eq_film_reg()

        flow_structure = self.flow_structure
       # H_L_L_S = self.H_L_L_S
       # H_L_total_an = self.calc_H_L_gas_core()     #Без учета холдапа в пленке в аннулар режиме
    #    H_L_total_an = calc_H_L_total_an(v_s_g, v_s_l)
       # H_L_dispers = self.calc_H_L_dispers()
        H_L_strat = self.H_L_strat #(self.dim_h_l)
        H_L_S_U = self.H_L_S_U
        H_L_T_B = self.calc_H_L_T_B_fact() #(self.v_s_g, self.v_s_l)
        d = [H_L_T_B,H_L_S_U,self.H_L_L_S]
        if flow_structure == 1:
            return 0# H_L_total_an
        elif flow_structure == 2:
            return H_L_strat
        elif flow_structure == 3:
            return H_L_strat
        elif flow_structure == 4:
            return self.H_L_dispers
        else:
            return d
    

    #Модуль расчета градиента давления
    #Расчет градиета давления
    def calc_pressure_drop_total(self):
        
        
        self.calc_auxiliary_pres_drop()
        #Intermittent Flow (Slug Regime)
        a1 = self.L_f/self.L_U
        b1 = self.rho_f*self.gravity*math.sin((math.pi*self.angle/180))
        c1 = (self.tau_w_film*self.S_f + self.tau_w_gas_bubble*self.S_g)/self.Area_pipe_m
        d1 = a1*(c1 + b1)
        self.pressure_drop_film_reg = -d1      #Pressure Drop in Film Region
        
        a2 = self.L_S/self.L_U
        b2 = self.tau_Slug*4/self.d_pipe
        c2 = self.rho_slug_PT*self.gravity*math.sin((math.pi*self.angle/180))
        d2 = a2*(b2 + c2)
        self.calc_pressure_drop_slug_reg = -d2     #Pressure Drop in Slug Region

        self.pressure_drop_S_U = self.calc_pressure_drop_slug_reg + self.pressure_drop_film_reg #Pressure Drop in Intermittent Flow
        
        #Dispersed Bubble Flow
        a3 = self.tau_no_slip*4/self.d_pipe
        b3 = self.rho_no_slip*self.gravity*math.sin((math.pi*self.angle/180))
        c3 = a3+b3
        self.pressure_drop_dispers = - c3  #Pressure Drop in Dispersed Bubble Flow
 
        #Stratified Flow
        a4 = self.tau_i_str*self.S_i_Str/self.Area_gas
        b4 = self.tau_w_g_str*self.S_g_Str/self.Area_gas
        c4 = self.rho_g_PT*self.gravity*math.sin((math.pi*self.angle/180))
        d4 = a4 + b4 + c4
        self.pressure_drop_gas_str = -d4

        a5 = -self.tau_i_str*self.S_i_Str/self.Area_liquid_Str
        b5 = self.tau_w_L_str*self.S_L_Str/self.Area_liquid_Str
        c5 = self.rho_l_PT*self.gravity*math.sin((math.pi*self.angle/180))
        d5 = a5 + b5 + c5
        self.pressure_drop_liquid_str = -d5
    
        self.pressure_drop_str = self.pressure_drop_liquid_str + self.pressure_drop_gas_str #Pressure Drop in Stratified Flow
    
        e = self.flow_structure
        if e == 2:
            return self.pressure_drop_str
        elif e == 3:
            return self.pressure_drop_str
        elif e == 4:
            return self.pressure_drop_dispers
        elif e == 5:
            return self.pressure_drop_S_U


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

grX = GradXiao()
grX.calc_auxiliary_init_vars()

print ("L_S = ", grX.L_S)
print ("L_S_max = ", grX.L_S_max)
#b= grX.calc_uniq_func_dim_height(0.01)
d= grX.calc_regime_liquid_holdup()

print ("Holdup =", d)
print (grX.calc_auxiliary_dimhl_vars())

print ("dP =", grX.calc_pressure_drop_total())
print ("dP_liq_str =", grX.pressure_drop_liquid_str)
print ("dp_gas_str =", grX.pressure_drop_gas_str)



"""
a = grX.combined_momentum_eq_film_reg()
b= grX.calc_regime_liquid_holdup()
print ("L_U = ", grX.calc_L_U(vg,vl,0.8))
print (grX.calc_pressure_drop_slug_reg(vg,vl, 0.8))
print (grX.calc_pressure_drop_film_reg(vg,vl, 0.8))
print (grX.combined_momentum_eq_film_reg())
print (grX.calc_H_L_T_B_fact(vg,vl))
print (grX.calc_pressure_drop_S_U(vg,vl, 0.8))
print (grX.pressure_drop(vg,vl, 0.8))


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
"""
#print (h_f_d)
"""
t=np.arange(0.2, 0.9, 0.01)
yy=[]
 
for xx in t :
    yy.append(grX.calc_uniq_func_dim_height(xx))
plt.plot(t,yy) 
plt.grid(True)
#plt.interactive(False)
plt.show()
"""
