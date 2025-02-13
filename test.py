import streamlit as st
import geemap

# Streamlit ilovasi sarlavhasi
st.title("Interaktiv Polygon Chizish va NDVI Vizualizatsiyasi")

# Geemap xaritasi
Map = geemap.Map(center=[40.7128, -74.0060], zoom=10)
