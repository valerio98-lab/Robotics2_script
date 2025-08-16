function check_rotation_matrix(R, checkRodriguez, theta, r)

        % Imposta la tolleranza predefinita se non viene fornita
        tol = 1e-10;

        if nargin < 2 || isempty(checkRodriguez)
            checkRodriguez = false;
        end

        if checkRodriguez && (nargin < 4 || isempty(theta) || isempty(r))
            error('Se checkRodriguez è true, è necessario fornire theta.');
        end
    
        % Controlla che R sia una matrice 3x3
        [m, n] = size(R);
        if m ~= 3 || n ~= 3
            warning('La matrice non è di dimensione 3x3.');
            isValid = false;
            return;
        end
    
        % Verifica l'ortogonalità: R'*R deve essere uguale alla matrice identità
        I = eye(3);
        orthoCheck = norm(R' * R - I, 'fro') < tol;
    
        % Verifica che il determinante di R sia 1
        detCheck = abs(det(R) - 1) < tol;
    
        % La matrice è una matrice di rotazione solo se entrambe le condizioni sono verificate
        isValid = orthoCheck && detCheck;

        if checkRodriguez
            traceVal = trace(R);
            expectedTrace = 1 + 2*cos(theta);
            rodCheck = abs(traceVal - expectedTrace) < tol;
            R_pos = double(RodriguezMatrix(r, theta));
            R_neg = double(RodriguezMatrix(-r, -theta));
            flag = isequal(R_pos, R_neg);
            if ~flag && rodCheck
                disp('Matrix checked has a valid trace but is not invariant')
            end
            % La matrice è valida solo se soddisfa anche questo controllo
            isValid = isValid && rodCheck && flag;
            
        end

        if isValid && checkRodriguez
            disp('is Valid and Rodriguez'); 
        elseif isValid
            disp('is Valid');
        else 
            disp('is not Valid. Watch Out!');
        end
    end
    
