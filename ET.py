def evapotranspiration(image):
        """ Evapotranspiratsiyani hisoblash funksiyasi """
        bands_needed = ['Rn', 'G', 'H']

        # Tasvirda kerakli bandlar borligini tekshiramiz
        missing_bands = [band for band in bands_needed if not image.bandNames().contains(band)]

        if missing_bands:
            raise ValueError(f"Quyidagi bandlar yetishmayapti: {missing_bands}")

        r_n = image.select('Rn')
        g = image.select('G')
        h = image.select('H')

        # ET = Rn - G - H
        et = r_n.subtract(g).subtract(h).rename('ET')

        return image.addBands(et)
