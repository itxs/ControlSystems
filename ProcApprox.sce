///////////////////////////////////////////////////////////////////////
/// Author: Khasanshin Timur Salavatovich, 2021
/// Licensed under GNU GPL v3
///
/// ProcApprox v1.0 - Library for system identification
/// 
/// ProcApprox is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// ProcApprox is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with ProcApprox. If not, see <https://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////

function y = FOPDT(x, k)
    for i=1:size(x,"*")
        if x(i) < k(2) then
            y(i) = 0
        else
            t = x(i)-k(2)
            a = k(1)
            b = k(3)
            y(i) = (1-exp(-t/b))/a
        end
    end
endfunction

function y = SOPDT(x, k)
    for i=1:size(x,"*")
        if x(i) < k(2) then
            y(i) = 0
        else
            t = x(i)-k(2)
            a = k(1)
            b = k(3)
            c = k(4)
            y(i) = (exp(-(b*t)/(2*a))*(-(a*b*sinh((sqrt(b^2-4*a*c)*t)/(2*a)))/(c*sqrt(b^2-4*a*c))-(a*cosh((sqrt(b^2-4*a*c)*t)/(2*a)))/c))/a+1/c
        end
    end
endfunction

function e=SOPDT_errorFunc(k, z)
    e=real(z(2,:) - SOPDT(z(1,:), k)) + imag(z(2,:) - SOPDT(z(1,:), k))
endfunction

function e=FOPDT_errorFunc(k, z)
    e=z(2,:) - FOPDT(z(1,:), k)
endfunction

function [k, e]=ProcessApproximate(x, y, process, startVals, weights)
    if ~exists('weights') then
        weights = ones(size(x,"*"),1)
    end
    if convstr(process) == "fopdt" then
        if ~exists('startVals') then
            startVals = ones(3,1)*5
        end
        [k, e] = datafit(FOPDT_errorFunc, [x;y], weights', startVals');
        plot2d(x', [FOPDT(x', k), y'])
    elseif convstr(process) == "sopdt" then
        if ~exists('startVals') then
            startVals = ones(4,1)*5
        end
        [k, e] = datafit(SOPDT_errorFunc, [x;y], weights', startVals');
        plot2d(x', [SOPDT(x', k), y'])
    else
        k = 0
        e = 0
        disp("Wrong input parameter.")
        return
    end
    disp("Coeffs: ", k)
    disp("Error sum: ", e)
endfunction
