function [output,e1,e2] = subDomainSpring(Force_pN,threshold)
%SUBDOMAINSPRING Summary of this function goes here
%   Detailed explanation goes here

k1 = 13;     % Firm Spring Constant representing Base Sub-Domain Elastisity
k2 = 1e-5;     % Softer Spring Constant representing the unfolding sections of the Sub-Domain

% foldedLength = 5;       % nm, Assumed Length
foldedLength = 0;
unfoldedLength = 40;    % nm, Assumed Length

% If Force = Constant * Elongation
% Elongation = Force / Constant

% I only expect the k1 spring to expand a couple of nm over the 0-30 pN
% range

e1 = 0;
e2 = 0;

a = 3e-1;
b = 1.5;
% e2 = nthroot((Force_pN/a),b);

if Force_pN >= threshold
    e2 = Force_pN / k2;
%     e2 = nthroot((Force_pN/a),b);
    if e2 >= unfoldedLength - foldedLength
        e2 = unfoldedLength - foldedLength;
    end
end

e1 = Force_pN / k1;

output = foldedLength + e1 + e2;

end

