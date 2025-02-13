def calc_aerodinamic_res(image, wind_image):
        """ AERODINAMIK QARSHILIKNI HISOBLAYMIZ """
        # Ma`lumootlarni boshlash
        lai = image.select('LAI')
        wind_speed = wind_image.select('u_component_of_wind_10m')  # Fixed variable name (wind_image instead of wind_speed)

        # Z_OM ni hisoblash
        z_om = ee.Image.constant(0.018).multiply(lai)

        # SHAMOL TEZLIGI VA DOIMIY
        wind_height = ee.Image.constant(10)
        k = ee.Image.constant(0.41)  # von Kármán constant

        # Friction velocity ni hisoblash
        u_f_v = k.multiply(wind_speed).divide(wind_height.divide(z_om).log())

        # Z lar qiymatini hisoblash
        z_2 = 2
        z_1 = 1

        # Aerodynamic resistance (R_AH) ni hisoblaash
        r_ah = image.expression(
            'log(z2 / z1) / (u * k)',  # Fixed log syntax for better clarity
            {
                'z2': z_2,
                'z1': z_1,
                'u': u_f_v,
                'k': k
            }
        ).rename('R_AH')

        return image.addBands(r_ah)
    
def different_temp(image, air_t):
        """ Surface temperature and air temperature lari o'rtasidagi farq"""
        s_t = image.select('T_S')  # Surface temperature
        a_t = air_t.select('temperature_2m')  # Air temperature 2m dagi
        d_t = s_t.subtract(a_t).rename('dT')  # Haroratlar farqi
        return image.addBands(d_t)  # yangi band qo'shish

def calc_havo_bosim(image, bosim, harorat):
        """ Havo zichligini hisoblash """
        p = bosim.select('mean_sea_level_pressure')  # Bosimni tanlash
        t = harorat.select('temperature_2m')  # Haroratni tanlash
        r_d = ee.Image.constant(287.5)  # Gaz konstantasi (havo uchun)

        # Havo zichligini hisoblash
        r = image.expression(
            'p / (R * t)',  # Bu yerda R va t o'zgaruvchilariga to'g'ri ishlov berildi
            {
                'p': p,
                't': t,
                'R': r_d
            }
        ).rename('Bosim')  # Havo zichligini nomlash

        return image.addBands(r)  # R bandini tasvirga qo'shish
    
def calc_heat_flux(image):
        """ H ni hisoblash funksiyasi """
        t = image.select('T_S')  # Yuzaning harorati
        r = image.select('R_AH')  # Aerodinamik qarshilik
        cp = ee.Image.constant(1004)  # Havo doimiysi (J/(kg·K))
        dt = image.select('dT')  # Harorat farqi (dT)

        h = image.expression(
            'cp * dT / rah',  # Heat flux formula
            {
                'cp': cp,
                'dT': dt,
                'rah': r
            }
        ).rename('H')  # Natijaviy bandni nomlash

        return image.addBands(h)  # H bandini tasvirga qo'shish
