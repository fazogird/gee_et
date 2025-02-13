def soil_heat_flux(image):
        """ G ni SOIL HEAT flux ni hisoblaymiz """
        r_n = image.select('Rn')
        t_s = image.select('T_S')
        a = image.select('ALBEDO')
        ndvi = image.select('NDVI')

        g = image.expression(
            '(Rn * Ts) / a * (0.0038 * a + 0.0074 * a.pow(2))) * (1 - 0.98 * ndvi.pow(4)',
            {
                'Rn': r_n,
                'Ts': t_s,
                'a': a,
                'ndvi': ndvi
            }
        ).rename('G')

        # band qo'shish
        return image.addBands(g)