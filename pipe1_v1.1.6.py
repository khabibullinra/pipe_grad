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
import unitconverter as un

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
un.g = 9.81  #ускорение свободного падения, м/с2
sigma_l = 0.03 #поверхностное натяжение на границе жидкости и газа, Н/м

"""


class Grad:
    def __init__(self):
        # Исходные данные для проведения расчета
        self.v_s_l = 0.2  # приведенная скорость жидкости, м3/с
        self.v_s_g = 1  # приведенная скорость газа, м3/с
        self.P_pipe = 40  # давление в трубе, МПа
        self.T_pipe = 350  # температура на устье скважины, К
        self.L_pipe = 2000  # длина трубы, м
        self.d_pipe = 0.8  # диаметр трубы, м
        self.rho_l_PT = 1000  # плотность жидкости при термобарических условиях, кг/м3
        self.rho_g_PT = 20  # плотность газа при термобарических условиях, кг/м3
        self.mu_l_PT = 0.001  # вязкость жидкости при термобарических условиях, Па*с
        self.mu_g_PT = 0.00002  # вязкость газа при термобарических условиях, Па*с
        self.angle = 0  # угол наклона трубы, градус
        self.rough = 0.0001  # абсолютная шероховатость трубы, м
        # self.un.g = un.g       # ускорение свободного падения, м/с2
        self.sigma_l = 0.03  # поверхностное натяжение на границе жидкости и газа, Н/м

    @property
    def v_s_mix(self):
        return self.v_s_g + self.v_s_l

    @property
    def n_re_liq(self):
        return self.rho_l_PT * self.d_pipe * self.v_s_l / self.mu_l_PT

    @property
    def n_re_gas(self):
        return self.rho_l_PT * self.d_pipe * self.v_s_g / self.mu_g_PT

    @property
    def n_re_mix_slug(self):
        return self.v_s_mix * self.rho_l_PT * self.d_pipe / self.mu_l_PT

    @property
    def area_pipe_m(self):
        return un.pi * self.d_pipe ** 2 / 4  # диаметр трубы


class GradXiao(Grad):
    
    # Расчет предварительных величин
    def __init__(self):
        """
        конструктор класса, задает исходные значения 
        """
<<<<<<< HEAD
        super().__init__()
        """ переменные вычисляемые в ходе расчетов """
        # self.area_pipe_m = 0
        self.l_s = 0
        self.l_s_max = 0
        self.N_Re_mix_slug_crit = 1000
        self.v_t_b = 0

=======
        # Исходные данные для проведения расчета
<<<<<<< HEAD
        self.v_s_l = 0.2          # приведенная скорость жидкости, м3/с
        self.v_s_g = 10           # приведенная скорость газа, м3/с
=======
        self.v_s_l = 2          # приведенная скорость жидкости, м3/с
        self.v_s_g = 20            # приведенная скорость газа, м3/с
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
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
        
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
    def calc_auxiliary_init_vars(self):  # расчет вспомогательных переменных
        # self.area_pipe_m =
        d_pipe_inch = self.d_pipe / 0.0254
        # расчет длины слага, зависит только от диаметра
        a = math.log(d_pipe_inch)
        b = -25.4 + 28.5 * a ** 0.1
        l_s_feet = math.exp(b)
        self.l_s = l_s_feet * 0.3048
        # расчет максимальной длины слага 
        b = math.log(l_s_feet)
        c = 0.5 * 3.08 + b
        l_s_feet_max = math.exp(c)
        self.l_s_max = l_s_feet_max * 0.3048
        # self.v_s_mix = self.v_s_g + self.v_s_l
        # найдем число Рейнольдса для жидкости для заданных условий
<<<<<<< HEAD
        # self.n_re_liq = self.rho_l_PT * self.d_pipe * self.v_s_l / self.mu_l_PT
        # self.n_re_gas = self.rho_g_PT * self.d_pipe * self.v_s_g / self.mu_g_PT
        c_o = 2.27 / (1 + (self.n_re_mix_slug / self.N_Re_mix_slug_crit) ** 2) \
              + 1.2 / (1 + (self.N_Re_mix_slug_crit / self.n_re_mix_slug) ** 2)

        aa = (un.g * self.d_pipe) ** 0.5
        v_d_tb = 0.54 * aa * math.cos(math.pi * self.angle / 180) \
                 + 0.35 * aa ** 0.5 * math.sin(math.pi * self.angle / 180)
        self.v_t_b = c_o * self.v_s_mix + v_d_tb
        self.H_L_L_S = 1 / (1 + (self.v_s_mix / 8.66) ** 1.39)
=======
        self.N_Re_l = self.rho_l_PT * self.d_pipe*self.v_s_l / self.mu_l_PT
        self.N_Re_g = self.rho_g_PT * self.d_pipe * self.v_s_g / self.mu_g_PT
        self.N_Re_mix_slug = self.v_s_mix * self.rho_l_PT * self.d_pipe / self.mu_l_PT
        self.N_Re_mix_slug_crit = 1000
        self.c_o = 2.27/(1+(self.N_Re_mix_slug / self.N_Re_mix_slug_crit)**2) + 1.2 / (1+(self.N_Re_mix_slug_crit / self.N_Re_mix_slug)**2)
<<<<<<< HEAD
        aa = (self.gravity*self.d_pipe)**0.5
        self.v_d_TB = 0.54 * aa * math.cos(math.pi * self.angle/180) + 0.35 * aa ** 0.5 * math.sin(math.pi * self.angle / 180)
=======
        a = (self.gravity*self.d_pipe)**0.5
        self.v_d_TB = 0.54 * a * math.cos(math.pi * self.angle/180) + 0.35 * a ** 0.5 * math.sin(math.pi * self.angle / 180)
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
        self.v_T_B = self.c_o * self.v_s_mix + self.v_d_TB
        self.H_L_L_S = 1/(1+(self.v_s_mix/8.66)**1.39)
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
        self.H_L_dispers = self.v_s_l / (self.v_s_l + self.v_s_g)
        c = (un.g * self.sigma_l * (self.rho_l_PT - self.rho_g_PT)
             / self.rho_l_PT ** 2) ** 0.25
        self.v_g_L_S = self.v_s_mix + 1.53 * c * self.H_L_L_S ** 0.1 \
                       * math.sin(math.pi * self.angle / 180)
        self.v_L_L_S = (self.v_s_l + self.v_s_g - self.v_g_L_S * (1 - self.H_L_L_S)) / self.H_L_L_S
        self.H_L_S_U = (self.v_t_b * self.H_L_L_S + self.v_g_L_S * (1 - self.H_L_L_S) - self.v_s_g) / self.v_t_b

        if self.n_re_liq < 2000:
            self.n_power_coeff = 1
<<<<<<< HEAD
            self.f_factor_l = 64 * self.n_re_liq ** (-1)
=======
<<<<<<< HEAD
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
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
        else:
            self.n_power_coeff = 0.2
            self.f_factor_l = 0.184 * self.n_re_liq ** (-0.2)

<<<<<<< HEAD
        if self.n_re_gas < 2000:
=======
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
=======
        else:
            self.n_power_coeff = 0.2       

        if self.N_Re_g<2000:
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
            self.m_power_coeff = 1
            self.f_factor_g = 64 * self.n_re_gas ** (-1)
        else:
            self.m_power_coeff = 0.2
            self.f_factor_g = 0.184 * self.n_re_gas ** (-0.2)

        self.pres_drop_s_l = - self.f_factor_l * self.rho_l_PT * self.v_s_l ** 2 / (2 * self.d_pipe)
        self.pres_drop_s_g = - self.f_factor_g * self.rho_g_PT * self.v_s_g ** 2 / (2 * self.d_pipe)
        self.x_large = ((-self.pres_drop_s_l) / (-self.pres_drop_s_g)) ** 0.5
        self.y_large = (self.rho_l_PT - self.rho_g_PT) * un.g * math.sin(math.pi * self.angle / 180) / (
        -self.pres_drop_s_g)
        self.f_large = (self.rho_g_PT / (self.rho_l_PT - self.rho_g_PT)) ** 0.5 * self.v_s_g / (
                                                                                               self.d_pipe * un.g * math.cos(
                                                                                                   math.pi * self.angle / 180)) ** 0.5
        self.k_large = self.n_re_liq ** 0.5 * self.f_large
        self.t_large = (-self.pres_drop_s_l / (
        (self.rho_l_PT - self.rho_g_PT) * un.g * math.cos(math.pi * self.angle / 180))) ** 0.5
        # Slug Liquid Holdup in Slug Region (fig 4.14)
        return

    def calc_auxiliary_hdf_vars(self):
        """
        Здесь расчитываются параметры, которые надо пересчитать после оценнки толщины пленки
        :return:
        """

        # Определение Liquid Holdup in Taylor Bubble
        # Вычисление холдапа жидкоти в Тейлор бабл регионе из предположения размера пленки, как для стратифайд режима
        self.H_L_T_B = 0.5 + (4 / math.pi) * (self.h_f_d - 0.5) * (self.h_f_d - self.h_f_d ** 2) ** 0.5 + (
                                                                                                              1 / math.pi) * math.asin(
            2 * self.h_f_d - 1)
        # Actual film velocity in Taylor Bubble Region
        self.v_L_T_B = self.v_t_b - (self.v_t_b - self.v_L_L_S) * self.H_L_L_S / self.H_L_T_B
        # Actual Taylor Bubble velocity in Taylor Bubble Region
        self.v_g_T_B = (self.v_s_mix - self.v_L_T_B * self.H_L_T_B) / (1 - self.H_L_T_B)
        # Геометрические параметры в для региона Тейлор бабл
        b = 2 * self.h_f_d - 1
        self.S_f = self.d_pipe * (math.pi - math.acos(b))
        self.S_g = math.pi * self.d_pipe - self.S_f
        b = 2 * self.h_f_d - 1
        self.S_i = self.d_pipe * (1 - b * b) ** 0.5
        self.Area_film = self.H_L_T_B * self.area_pipe_m
        self.Area_gas = (1 - self.H_L_T_B) * self.area_pipe_m
        self.d_hydr_film = 4 * self.Area_film / self.S_f  # расчет гидравлического диаметра для пленки
        self.d_hydr_gas_bubble = 4 * self.Area_gas / (
        self.S_g + self.S_i)  # расчет гидравлического диаметра для пузырька Тейлора
        self.N_Re_film = self.rho_l_PT * self.d_hydr_film * math.fabs(
            self.v_L_T_B) / self.mu_l_PT  # Расчет числа Рейнольдса для пленки
        self.N_Re_gas_bubble = self.rho_g_PT * self.d_hydr_gas_bubble * math.fabs(
            self.v_g_T_B) / self.mu_g_PT  # Расчет числа Рейнольдса для пузырька Тейлора
        if self.N_Re_film < 2000:  # Расчет гидравлического коэффициента потерь на трение для пленки
            self.f_factor_film = 64 * self.N_Re_film ** (-1)  # режим ламинарный
        else:
            self.f_factor_film = 0.184 * self.N_Re_film ** (-0.2)  # режим турбулентный
        # Расчет гидравлического коэффициента потерь на трение для пузырька Тейлора
        if self.N_Re_gas_bubble < 2000:
            self.f_factor_gas_bubble = 64 * self.N_Re_gas_bubble ** (-1)  # режим ламинарный
        else:
            self.f_factor_gas_bubble = 0.184 * self.N_Re_gas_bubble ** (-0.2)  # режим турбулентный
        # Shear Stress
        # Liquid Wall Shear Stress for film in Taylor Bubble Region
        self.tau_w_film = self.f_factor_film * self.rho_l_PT * self.v_L_T_B * math.fabs(self.v_L_T_B) / 8
        #  Gas Wall Shear Stress for gas bubble in Taylor Bubble Region
        self.tau_w_gas_bubble = self.f_factor_gas_bubble * self.rho_g_PT * self.v_g_T_B * math.fabs(self.v_g_T_B) / 8
        # Interfacial Shear Stress
        f_factor_i = 0.0568
        self.tau_i = f_factor_i * self.rho_g_PT * (self.v_g_T_B - self.v_L_T_B) * math.fabs(
            self.v_g_T_B - self.v_L_T_B) / 8
        a = self.l_s * (self.v_L_L_S * self.H_L_L_S - self.v_L_T_B * self.H_L_T_B)
        b = self.v_s_l - self.v_L_T_B * self.H_L_T_B
        self.L_U = a / b
        self.L_f = self.L_U - self.l_s
        self.Freq_S = self.v_t_b / self.L_U
        # Empirical correlation for Slug Frequency by Zabaras (2000)
        v_s_l_feet = self.v_s_l / 0.3048
        un.g_feet = 32.2
        d_pipe_feet = self.d_pipe / 0.3048
        v_mix_feet = self.v_s_mix / 0.3048
        a = (v_s_l_feet / (un.g_feet * d_pipe_feet)) ** 1.2
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
        s = 0.01
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
<<<<<<< HEAD
=======
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
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
        elif trans_A >= 1 and trans_B > 0.35:
            nonstratified_regime = bubble_slug_regime
            # 1-Annular; 2-Stratified-Wavy; 3-Stratified-Smooth; 4-Dispersed-Bubble; 5-Intermittent
            #    trans_A = self.calc_trans_A(dim_h_l)
        if trans_A < 1:
            self.flow_structure = stratified_regime
        else:
            self.flow_structure = nonstratified_regime
<<<<<<< HEAD

        # Stratified flow
        self.H_L_strat = 0.5 + (4 / math.pi) * (self.dim_h_l - 0.5) * (self.dim_h_l - self.dim_h_l ** 2) ** 0.5 + (1
                                                                                                                   / math.pi) * math.asin(
            2 * self.dim_h_l - 1)
        return self.flow_structure

    def calc_auxiliary_pres_drop(self):
        # Slug Flow
        H_L_T_B = self.calc_H_L_T_B_fact()
        self.rho_f = self.rho_l_PT * H_L_T_B + self.rho_g_PT * (1 - H_L_T_B)
        self.mu_slug_PT = self.mu_l_PT * self.H_L_L_S + self.mu_g_PT * (1 - self.H_L_L_S)  # Вязкость слага
        self.rho_slug_PT = self.rho_l_PT * self.H_L_L_S + self.rho_g_PT * (1 - self.H_L_L_S)  # Плотность слага
        self.N_Re_Slug = self.rho_slug_PT * self.d_pipe * self.v_s_mix / self.mu_slug_PT  # Расчет числа Рейнольдса для Slug
        # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.N_Re_Slug < 2000:
            # режим ламинарный
            self.f_factor_Slug = 64 * self.N_Re_Slug ** (-1)
        else:
            # режим турбулентный
            self.f_factor_Slug = 0.184 * self.N_Re_Slug ** (-0.2)

        self.tau_Slug = self.f_factor_Slug * self.rho_slug_PT * self.v_s_mix * self.v_s_mix / 8  # Shear Stress in Slug Region

        # Dispersed Bubble Flow
        self.rho_no_slip = self.rho_l_PT * self.H_L_dispers + self.rho_g_PT * (1 - self.H_L_dispers)
        self.mu_no_slip = self.H_L_dispers * self.mu_l_PT + (1 - self.H_L_dispers) * self.mu_g_PT
        self.Re_no_slip = self.v_s_mix * self.d_pipe * self.rho_no_slip / self.mu_no_slip
        # Расчет гидравлического коэффициента потерь на трение для Slug
        if self.Re_no_slip < 2000:
            # режим ламинарный
            self.f_factor_no_slip = 64 * self.Re_no_slip ** (-1)
        else:
            # режим турбулентный
            self.f_factor_no_slip = 0.184 * self.Re_no_slip ** (-0.2)

        self.tau_no_slip = self.f_factor_no_slip * self.rho_no_slip * self.v_s_mix * self.v_s_mix / 8  # Shear Stress in Dispersed Region

        # Stratified Flow (Smooth and Wavy)
        h_l = self.calc_dim_h_l_fact(self.dim_h_l)
        b = 2 * h_l - 1
        self.S_L_Str = self.d_pipe * (math.pi - math.acos(b))
        self.S_g_Str = math.pi * self.d_pipe - self.S_L_Str
        self.S_i_Str = self.d_pipe * (1 - b * b) ** 0.5
        self.Area_liquid_Str = self.H_L_strat * self.area_pipe_m
        self.Area_gas_Str = (1 - self.H_L_strat) * self.area_pipe_m
        self.d_hydr_liq = 4 * self.Area_liquid_Str / self.S_L_Str  # ПРОВЕРИТЬ! расчет гидравлического диаметра для слоя жидкости
        self.d_hydr_gas = 4 * self.Area_gas / (
        self.S_g_Str + self.S_i_Str)  # ПРОВЕРИТЬ #расчет гидравлического диаметра для слоя газа
        self.v_L_Str = self.v_s_l / self.H_L_strat  # Actual liquid velocity
        self.v_g_Str = self.v_s_g / (1 - self.H_L_strat)  ##Actual gas velocity
        self.N_Re_gas_Str = self.v_g_Str * self.d_hydr_gas * self.rho_g_PT / self.mu_g_PT
        self.N_Re_liquid_Str = self.v_L_Str * self.d_hydr_liq * self.rho_l_PT / self.mu_l_PT
        # Расчет гидравлического коэффициента потерь на трение для жидкости
        if self.N_Re_liquid_Str < 2000:
            # режим ламинарный
            self.f_factor_l_str = 64 * self.N_Re_liquid_Str ** (-1)
        else:
            # режим турбулентный
            self.f_factor_l_str = 0.184 * self.N_Re_liquid_Str ** (-0.2)
        # Расчет гидравлического коэффициента потерь на трение для газа
        if self.N_Re_gas_Str < 2000:
            # режим ламинарный
            self.f_factor_g_str = 64 * self.N_Re_gas_Str ** (-1)
        else:
            # режим турбулентный
            self.f_factor_g_str = 0.184 * self.N_Re_gas_Str ** (-0.2)
        # Расчет гидравлического коэффициента потерь на трение между жидкостью и газом
        aa = self.flow_structure
        cc = 0.0568
        if aa == 2:
            self.f_factor_i_str = cc
        elif aa == 3:
            self.f_factor_i_str = self.f_factor_g_str
        else:
            self.f_factor_i_str = 0

        self.tau_w_L_str = self.f_factor_l_str * self.rho_l_PT * self.v_L_Str * self.v_L_Str / 8
        self.tau_w_g_str = self.f_factor_g_str * self.rho_g_PT * self.v_g_Str * self.v_g_Str / 8
        self.tau_i_str = 0.125 * self.f_factor_i_str * self.rho_g_PT * (self.v_g_Str - self.v_L_Str) * (
        self.v_g_Str - self.v_L_Str)

        return
=======
<<<<<<< HEAD

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


=======

        # Stratified flow
        self.H_L_strat = 0.5 + (4 / math.pi) * (self.dim_h_l - 0.5) * (self.dim_h_l - self.dim_h_l ** 2) ** 0.5 + (1
                                                                   / math.pi) * math.asin(2 * self.dim_h_l - 1)

        return  self.flow_structure


>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4

    def calc_uniq_func_dim_height(self, dim_h_l):
        def calc_dim_A_l(dim_h_l):
            b = 2 * dim_h_l - 1
            a = math.pi - math.acos(b)
            c = 1 - b * b
            self.dim_A_l = 0.25 * (a + b * c ** 0.5)
            return self.dim_A_l

        def calc_dim_A_g(dim_h_l):
            b = 2 * dim_h_l - 1
            a = math.acos(b)
            c = 1 - b ** 2
            self.dim_A_g = 0.25 * (a - b * c ** 0.5)
            return self.dim_A_g

        def calc_dim_S_l(dim_h_l):
            self.dim_S_l = math.pi - math.acos(2 * dim_h_l - 1)
            return self.dim_S_l

        def calc_dim_S_i(dim_h_l):
            self.dim_S_i = (1 - (2 * dim_h_l - 1) ** 2) ** 0.5
            return self.dim_S_i

        def calc_dim_S_g(dim_h_l):
            self.dim_S_g = math.acos(2 * dim_h_l - 1)
            return self.dim_S_g

        def calc_dim_v_l(dim_h_l):
            dim_A_g = 0.25 * (math.acos(2 * dim_h_l - 1) - (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            dim_A_l = 0.25 * (
            math.pi - math.acos(2 * dim_h_l - 1) + (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            self.dim_v_l = (dim_A_g / dim_A_l + 1)
            return self.dim_v_l

        def calc_dim_v_g(dim_h_l):
            dim_A_g = 0.25 * (math.acos(2 * dim_h_l - 1) - (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            dim_A_l = 0.25 * (
            math.pi - math.acos(2 * dim_h_l - 1) + (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            self.dim_v_g = 1 + dim_A_l / dim_A_g
            return self.dim_v_g

        def calc_dim_d_l(dim_h_l):
            dim_A_l = 0.25 * (
            math.pi - math.acos(2 * dim_h_l - 1) + (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            dim_S_l = math.pi - math.acos(2 * dim_h_l - 1)
            self.dim_d_l = 4 * dim_A_l / dim_S_l
            return self.dim_d_l

        def calc_dim_d_g(dim_h_l):
            dim_A_g = 0.25 * (math.acos(2 * dim_h_l - 1) - (2 * dim_h_l - 1) * (1 - (2 * dim_h_l - 1) ** 2) ** 0.5)
            dim_S_g = math.acos(2 * dim_h_l - 1)
            dim_S_i = (1 - (2 * dim_h_l - 1) ** 2) ** 0.5
            self.dim_d_g = 4 * dim_A_g / (dim_S_i + dim_S_g)
            return self.dim_d_g
<<<<<<< HEAD

        a1 = self.x_large
        a2 = self.y_large

        dim_h_l = self.calc_dim_h_l_fact(dim_h_l)
        b1 = calc_dim_A_l(dim_h_l)
        b2 = calc_dim_A_g(dim_h_l)
        c1 = calc_dim_S_l(dim_h_l)
        c2 = calc_dim_S_i(dim_h_l)
        c3 = calc_dim_S_g(dim_h_l)
        d1 = calc_dim_v_l(dim_h_l)
        d2 = calc_dim_v_g(dim_h_l)
        e1 = calc_dim_d_l(dim_h_l)
        e2 = calc_dim_d_g(dim_h_l)
        f1 = self.n_power_coeff
        f2 = self.m_power_coeff
        a11 = a1 * a1
        b11 = d1 * e1
        c11 = b11 ** (-f1)
        d11 = d1 * d1
        e11 = c1 / b1
        f11 = d2 * e2
        g11 = f11 ** (-f2)
        h11 = d2 * d2
        i11 = c3 / b2 + c2 / b1 + c2 / b2
        res = a11 * c11 * d11 * e11 - g11 * h11 * i11 + 4 * a2
        self.dim_h_l = dim_h_l
        return res

    def combined_momentum_eq_film_reg(self):  # Combined momentum equation in film region
        self.h_f_d = 0.001
        res = 0.1
        while (res > 0):
            #        h_f_d=h_f_d+0.0001
            self.calc_auxiliary_hfd_vars()
            tau_w_film = self.tau_w_film
=======
        
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






    #Dispersed Bubble Flow

    
    #Intermittent Flow
    #    def calc_v_s(self,v_s_g, v_s_l):         #Slug Velosity - средняя скорость гомогенной смеси в Slug Region (fig 4.14)
    #        return self.v_s_mix
    
    
   
    



            

>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21

    def combined_momentum_eq_film_reg(self):       #Combined momentum equation in film region
        self.h_f_d = 0.001
        res = 0.1
        while (res>0):
#        h_f_d=h_f_d+0.0001
<<<<<<< HEAD
            self.calc_auxiliary_hfd_vars()
=======
            self.calc_auxiliary_hdf_vars()
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
            tau_w_film= self.tau_w_film
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
            tau_w_gas_bubble = self.tau_w_gas_bubble
            tau_i = self.tau_i
            S_f = self.S_f
            S_g = self.S_g
            S_i = self.S_i
            Area_film = self.Area_film
            Area_gas = self.Area_gas
            a1 = tau_w_gas_bubble * S_g / Area_gas
            b1 = tau_w_film * S_f / Area_film
            c11 = tau_i * S_i
            c12 = 1 / Area_gas + 1 / Area_film
            c1 = c11 * c12
            d1 = (self.rho_l_PT - self.rho_g_PT) * un.g * math.sin(math.pi * self.angle / 180)
            # return b1
            res = a1 - b1 + c1 - d1
            self.h_f_d = self.h_f_d + 0.00001
        return res

    def calc_H_L_T_B_fact(self):
        if self.h_f_d >= 1:
            return 0.999
        elif self.h_f_d <= 0:
            return 0.001
        else:
            return self.h_f_d

    # Определение Liquid Holdup в замисимости от структуры потока
    def calc_regime_liquid_holdup(self):

        self.dim_h_l = fsolve(self.calc_uniq_func_dim_height, 0.001)

        self.calc_auxiliary_dimhl_vars()

        self.combined_momentum_eq_film_reg()

        flow_structure = self.flow_structure
        # H_L_L_S = self.H_L_L_S
        # H_L_total_an = self.calc_H_L_gas_core()     #Без учета холдапа в пленке в аннулар режиме
        #    H_L_total_an = calc_H_L_total_an(v_s_g, v_s_l)
        # H_L_dispers = self.calc_H_L_dispers()
        H_L_strat = self.H_L_strat  # (self.dim_h_l)
        H_L_S_U = self.H_L_S_U
        H_L_T_B = self.calc_H_L_T_B_fact()  # (self.v_s_g, self.v_s_l)
        d = [H_L_T_B, H_L_S_U, self.H_L_L_S]
        if flow_structure == 1:
            return 0  # H_L_total_an
        elif flow_structure == 2:
            return H_L_strat
        elif flow_structure == 3:
            return H_L_strat
        elif flow_structure == 4:
            return self.H_L_dispers
        else:
            return d

<<<<<<< HEAD
    # Модуль расчета градиента давления
    # Расчет градиета давления
    def calc_pressure_drop_total(self):

        self.calc_auxiliary_pres_drop()
        # Intermittent Flow (Slug Regime)
        a1 = self.L_f / self.L_U
        b1 = self.rho_f * un.g * math.sin((math.pi * self.angle / 180))
        c1 = (self.tau_w_film * self.S_f + self.tau_w_gas_bubble * self.S_g) / self.area_pipe_m
        d1 = a1 * (c1 + b1)
        self.pressure_drop_film_reg = -d1  # Pressure Drop in Film Region

        a2 = self.l_s / self.L_U
        b2 = self.tau_Slug * 4 / self.d_pipe
        c2 = self.rho_slug_PT * un.g * math.sin((math.pi * self.angle / 180))
        d2 = a2 * (b2 + c2)
        self.calc_pressure_drop_slug_reg = -d2  # Pressure Drop in Slug Region

        self.pressure_drop_S_U = self.calc_pressure_drop_slug_reg + self.pressure_drop_film_reg  # Pressure Drop in Intermittent Flow

        # Dispersed Bubble Flow
        a3 = self.tau_no_slip * 4 / self.d_pipe
        b3 = self.rho_no_slip * un.g * math.sin((math.pi * self.angle / 180))
        c3 = a3 + b3
        self.pressure_drop_dispers = - c3  # Pressure Drop in Dispersed Bubble Flow

        # Stratified Flow
        a4 = self.tau_i_str * self.S_i_Str / self.Area_gas
        b4 = self.tau_w_g_str * self.S_g_Str / self.Area_gas
        c4 = self.rho_g_PT * un.g * math.sin((math.pi * self.angle / 180))
        d4 = a4 + b4 + c4
        self.pressure_drop_gas_str = -d4

        a5 = -self.tau_i_str * self.S_i_Str / self.Area_liquid_Str
        b5 = self.tau_w_L_str * self.S_L_Str / self.Area_liquid_Str
        c5 = self.rho_l_PT * un.g * math.sin((math.pi * self.angle / 180))
        d5 = a5 + b5 + c5
        self.pressure_drop_liquid_str = -d5

        self.pressure_drop_str = self.pressure_drop_liquid_str + self.pressure_drop_gas_str  # Pressure Drop in Stratified Flow

=======
<<<<<<< HEAD
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
    
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
        e = self.flow_structure
        if e == 2:
            return self.pressure_drop_str
        elif e == 3:
            return self.pressure_drop_str
        elif e == 4:
            return self.pressure_drop_dispers
        elif e == 5:
            return self.pressure_drop_S_U

<<<<<<< HEAD
=======

=======

    

    #Модуль расчета градиента давления
    
    #Intermittent Flow (Slug Regime)
    def calc_pressure_drop_film_reg(self,v_s_g, v_s_l, d_pipe):      #Pressure Drop in Film Region
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        L_f = self.calc_L_f(v_s_g, v_s_l, d_pipe)
        tau_w_film= self.calc_tau_w_film(v_s_g, v_s_l)
        tau_w_gas_bubble = self.calc_tau_w_gas_bubble(v_s_g, v_s_l)
        S_f = self.calc_S_f()
        S_g = self.calc_S_g()
        Area_pipe = self.Area_pipe_m #calc_Area_pipe(d_pipe)
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
 #       v_s_mix = self.calc_v_s_mix(v_s_l, v_s_g)
        rho_slug = self.calc_rho_slug_PT(v_s_g, v_s_l)
        mu_slug_PT = self.calc_mu_slug_PT(self.mu_l_PT, self.mu_g_PT, v_s_g, v_s_l)
        return rho_slug*self.d_pipe*self.v_s_mix/mu_slug_PT
    
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
        v_s_m = self.v_s_mix
        return f_factor_Slug*rho_slug*v_s_m*v_s_m/8
    
    def calc_pressure_drop_slug_reg(self,v_s_g, v_s_l, d_pipe):      #Pressure Drop in Slug Region
        tau_Slug = self.calc_tau_Slug(v_s_g, v_s_l)
        L_U = self.calc_L_U(v_s_g, v_s_l, d_pipe)
        L_S = self.L_S
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
       # H_L_dispers = self.calc_H_L_dispers(v_s_g, v_s_l)
        return self.rho_l_PT * self.H_L_dispers + self.rho_g_PT*(1 - self.H_L_dispers)
    
    def calc_mu_no_slip(self,mu_l_PT, mu_g_PT):
        #H_L_dispers = self.calc_H_L_dispers(self.v_s_g, self.v_s_l) 
        return self.H_L_dispers*mu_l_PT + (1-self.H_L_dispers)*mu_g_PT
    
    def calc_Re_no_slip(self,v_s_g, v_s_l, d_pipe):
        rho_no_slip = self.calc_rho_no_slip(v_s_g, v_s_l)
        mu_no_slip = self.calc_mu_no_slip(self.mu_l_PT, self.mu_g_PT)
        v_s_m = self.v_s_mix
        return v_s_m*d_pipe*rho_no_slip/mu_no_slip
    
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
        v_s_m = self.v_s_mix
        return f_factor_no_slip*rho_no_slip*v_s_m*v_s_m/8
    
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
        Area_pipe = self.Area_pipe_m #calc_Area_pipe(self.d_pipe)
        H_L_strat = self.calc_H_L_strat(dim_h_l)
        return H_L_strat*Area_pipe
    
    def calc_Area_gas_Str(self,dim_h_l):
        Area_pipe = self.Area_pipe_m # calc_Area_pipe(self.d_pipe)
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
    
    def calc_pressure_drop_liquid_str(self):
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
    
    def calc_pressure_drop_str(self):
        a = self.calc_pressure_drop_liquid_str(v_s_g, v_s_l, d_pipe)
        b = self.calc_pressure_drop_gas_str(v_s_g, v_s_l, d_pipe)
        return a+b
    
    #Расчет градиета давления
    def pressure_drop(self):
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
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
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

<<<<<<< HEAD

grX = GradXiao()
grX.calc_auxiliary_init_vars()
=======
grX = GradXiao()
<<<<<<< HEAD
grX.calc_auxiliary_init_vars()
=======
grX.d_pipe = 0.8
grX.calc_auxiliary_init_vars()
dP = 0.8
vl = 1
vg = 55

grX.v_s_g = vg
grX.v_s_l = vl 
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21

print ("L_S = ", grX.L_S)
print ("L_S_max = ", grX.L_S_max)
#b= grX.calc_uniq_func_dim_height(0.01)
d= grX.calc_regime_liquid_holdup()
<<<<<<< HEAD

print ("Holdup =", d)
print (grX.calc_auxiliary_dimhl_vars())

print ("dP =", grX.calc_pressure_drop_total())
print ("dP_liq_str =", grX.pressure_drop_liquid_str)
print ("dp_gas_str =", grX.pressure_drop_gas_str)

=======
print(d)
>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21

>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4

print("l_s = ", grX.l_s)
print("l_s_max = ", grX.l_s_max)
# b= grX.calc_uniq_func_dim_height(0.01)
d = grX.calc_regime_liquid_holdup()

print("Holdup =", d)
print(grX.calc_auxiliary_dimhl_vars())

print("dP =", grX.calc_pressure_drop_total())
print("dP_liq_str =", grX.pressure_drop_liquid_str)
print("dp_gas_str =", grX.pressure_drop_gas_str)

vg_arr = np.arange(0.01,10.1,1)
vl_arr = np.arange(0.01,10.1,1)
xr = np.zeros([vl_arr.size, vg_arr.size])


for vg in vg_arr:
    grX.v_s_g = vg
    i=0
    for vl in vl_arr:
        grX.v_s_l = vl
        j=0
        xr1 = grX.calc_auxiliary_dimhl_vars()
        #grX.calc_auxiliary_dimhl_vars()
        #        print(grX.calc_regime_liquid_holdup())
        grX.calc_regime_liquid_holdup()
        print(grX.calc_regime_liquid_holdup())
        xr[i,j]=xr1
        j=j+1
    i=i+1
    
plt.pcolormesh(vl_arr, vg_arr, xr, cmap='RdBu', shading='flat')
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
# print (h_f_d)
"""
t=np.arange(0.2, 0.9, 0.01)
yy=[]
 
for xx in t :
    yy.append(grX.calc_uniq_func_dim_height(xx))
plt.plot(t,yy) 
plt.grid(True)
#plt.interactive(False)
plt.show()
<<<<<<< HEAD
"""
=======
<<<<<<< HEAD
"""
=======

>>>>>>> e325e29f960666c907b08a6fdf02ae997b902b21
>>>>>>> f0c91141a541fa6ac7f03d614fdb2eaa133038f4
