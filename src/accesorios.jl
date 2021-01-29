module accesorios

using Unitful

"""
    entrada_recta(m="Rennels")

Devuelve el coeficiente de pérdida para una entrada afilada a una 
tubería. Hay seis fuentes disponibles; cuatro de ellas recomiendan
K = 0.5, el más reciente método "Rennels" recomienda K = 0.57, y el 
método "Miller" recomienda de K = 0.51, leído en un gráfico.

![sdsdsd](assets/entrada_recta.png#thumbnail)

---

**Los campos de entrada son:**
- `m :: String, opcional ` 

Es el metodo a usar, uno de Rennels, Swamee, Blevins, Crane, Idelchik y Miller

---

**Salida:**
- `K :: Float ` 
Coeficiente de perdida `K`

"""
function entrada_recta(m="Rennels")
    if m==="Rennels"
        0.57
    elseif m in ("Swamee", "Blevins", "Crane", "Idelchik")
        0.5
    elseif m==="Miller"
        0.5092676683721356
    else
        error("Los métodos validos son: Rennels, Swamee, Blevins, Crane, Idelchik y Miller")
    end
end    

end

