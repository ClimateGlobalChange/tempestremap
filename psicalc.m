function psi = psicalc(dLonT, dLatT)

dLatC = (180/pi)*0.6;
dLonC = (180/pi)*0;

dSinC = sind(dLatC);
dCosC = cosd(dLatC);
dSinT = sind(dLatT);
dCosT = cosd(dLatT);

dTrm = dCosT .* cosd(dLonT - dLonC);
dX = dSinC .* dTrm - dCosC .* dSinT;
dY = dCosT .* sind(dLonT - dLonC);
dZ = dSinC .* dSinT + dCosC .* dTrm;

dLon = atan2(dY, dX);
for i = 1:length(dLon)
    if dLon(i) < 0
        dLon(i) = dLon(i) + 2*pi;
    end
end

dLat = asin(dZ);

dR0 = 3; dD = 5; dT = 6;

dRho = dR0 * cos(dLat);
dVt = (3*sqrt(3)/2) .* ((sech(dRho)).^2) .* tanh(dRho);
dOmega = zeros(size(dVt));

for i = 1:length(dRho)
    if dRho(i) == 0
        dOmega(i) = 0;
    else
        dOmega(i) = dVt(i)/dRho(i);
    end
end

psi = 1 - tanh((dRho/dD).*sin(dLon - dOmega * dT));

end
