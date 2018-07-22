function mgf_output = mgf(dist, phi, m, EbN0)
    if dist == 'naka'
        mgf_output = (1+EbN0./(m.*(sin(phi)).^2)).^(-m);
    else
        mgf_output = (1+EbN0./(sin(phi).^2)).^(-1);
    end
end