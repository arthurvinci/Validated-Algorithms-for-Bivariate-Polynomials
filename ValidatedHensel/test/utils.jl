function validate_error(F::Generic.Poly{Generic.AbsSeries{T}}, bound::Float64) where {T}
    for i in 0:(F.length - 1)
        sub_poly = coeff(F, i)
        for j in 0:(sub_poly.length - 1)
            if abs(Float64(coeff(sub_poly, j))) > bound
                return false
            end
        end
    end
    return true
end