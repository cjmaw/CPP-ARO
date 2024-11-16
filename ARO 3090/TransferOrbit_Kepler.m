rLEO = 550+6378; % km
rGEO = 36000+6378; % km
aXfer = (rLEO+rGEO)/2; % km

eXfer = (rGEO-rLEO)/(rGEO+rLEO);
[E,theta,M] = AngleCalculator(aXfer,eXfer);
r_theta=(aXfer.*(1-eXfer.^2))./(1+(eXfer.*cos(theta)));

endIndex_kepler = round(length(r_theta)/2); % To only plot the first half of the Xfer orbit
kepler_x = r_theta.*cos(theta);
kepler_x = kepler_x(1:endIndex_kepler);
kepler_y = r_theta.*sin(theta);
kepler_y = kepler_y(1:endIndex_kepler);

% ------------------------------------------------------------------------

function [EArray, thetaArray, MArray] = AngleCalculator(a,e)

    tolerance = 0.0001; 
    mu = 398600.435; % for Earth
    period = 2*pi*sqrt(a^3/mu); % seconds

    for t = 1:round(period)
        % Calculate Mean anomaly
        MArray(t) = sqrt(mu/(a^3))*t;
    
        % Newton's Method to Calculate Eccentric anomaly
        E(1) = MArray(t);
        f = @(E) E-e.*sin(E)-MArray(t);
        fprime = @(E) 1-e.*cos(E);
        
        for i=1:10
            E(i+1) = E(i) - (f(E(i))/fprime(E(i)));
            error(i) = abs(E(i)-E(i+1));
            if error(i)<tolerance
                break
            end
        end
        numIterations = i;
        EArray(t) = E(i);

        % Calculate True anomaly
        root = sqrt( (1+e)/(1-e) );
        thetaArray(t) = 2*atan(root*tan(EArray(t)/2));

    end

end
