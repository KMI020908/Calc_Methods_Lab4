#include "methods.h"

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        for (std::size_t i = k + 1; i < rows; i++){   // Частичный выбор главного элемента
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк       
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG gaussMethodFull(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    std::vector<std::size_t> switchCols; // Вектор, хранящий индексы перемещенных столбцов исходной матрицы
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        std::size_t mainСol = k; // Столбец главного элемента
        // Полный выбор
        // Поиск главного элмента в k-ом столбце  
        for (std::size_t i = k + 1; i < rows; i++){ 
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        // Поиск главного элмента в k-ой строке 
        for (std::size_t j = k + 1; j < cols; j++){ 
            if (std::abs(lCoefs[k][j]) > std::abs(mainValue)){
                mainValue = lCoefs[k][j]; 
                mainСol = j;
            }
        }
        //Замена столбцов
        if (mainСol != k){ 
            Type temp;
            for (std::size_t i = 0; i < rows; i++){
                temp = lCoefs[i][k];
                lCoefs[i][k] = lCoefs[i][mainСol];
                lCoefs[i][mainСol] = temp;
            }
            switchCols.push_back(k);
            switchCols.push_back(mainСol);
        }

        // Прямой ход Гаусса 
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    // Обратная перестановка
    for (int i = switchCols.size() - 2; i >= 0; i -= 2){
        Type temp = solution[switchCols[i]];
        solution[switchCols[i]] = solution[switchCols[i + 1]];
        solution[switchCols[i + 1]] = temp;

    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    if (rows != cols){
        solution.resize(0);
        return NO_SOLUTION;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/std::sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/std::sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c * lCoefs[k][j] + s * lCoefs[i][j];
                    lCoefs[i][j] = -s * temp + c * lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0;
                }
                Type temp = rCoefs[k];
                rCoefs[k] = c*rCoefs[k] + s*rCoefs[i];
                rCoefs[i] = -s*temp + c*rCoefs[i];
            }
        }
    }
    if (std::abs(lCoefs[rows - 1][rows - 1]) < accuracy){  // detA = 0
        return NO_SOLUTION;
    }
    
     // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }   
    return HAS_SOLUTION;
}

template<typename Type>
QUADRATIC_FLAG findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                Type s = lCoefs[i][k]/std::sqrt(lCoefs[k][k] * lCoefs[k][k] + lCoefs[i][k] * lCoefs[i][k]);
                for (std::size_t j = 0; j < cols; j++){
                    Type temp = Q[k][j];
                    Q[k][j] = c * Q[k][j] + s * Q[i][j];
                    Q[i][j] = -s * temp + c * Q[i][j];
                    if (std::abs(Q[i][j]) < accuracy)
                        Q[i][j] = 0.0;
                }
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c * lCoefs[k][j] + s * lCoefs[i][j];
                    lCoefs[i][j] = -s * temp + c * lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0.0;
                }
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findQMatrix3Diag(std::vector<std::vector<Type>> &matrix, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols);
        for (std::size_t j = 0; j < cols; j++){
            Q[i][j] = 0.0;   
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1.0;
    }
    // k < rows - 2
    for (std::size_t k = 0; k < rows - 2; k++){
        if (std::abs(matrix[k + 1][k]) >= accuracy){
            Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
            for (std::size_t j = 0; j < cols; j++){
                Type temp = Q[k][j];
                Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
                Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
                if (std::abs(Q[k + 1][j]) < accuracy)
                    Q[k + 1][j] = 0.0;
            }
            for (std::size_t j = k; j < k + 3; j++){
                Type temp = matrix[k][j];
                matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
                matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
                if (std::abs(matrix[k + 1][j]) < accuracy)
                    matrix[k + 1][j] = 0.0;
            }
        }
    }
    // k = rows - 2
    std::size_t k = rows - 2;
    if (std::abs(matrix[k + 1][k]) >= accuracy){
        Type c = matrix[k][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        Type s = matrix[k + 1][k]/std::sqrt(matrix[k][k] * matrix[k][k] + matrix[k + 1][k] * matrix[k + 1][k]);
        for (std::size_t j = 0; j < cols; j++){
            Type temp = Q[k][j];
            Q[k][j] = c * Q[k][j] + s * Q[k + 1][j];
            Q[k + 1][j] = -s * temp + c * Q[k + 1][j];
            if (std::abs(Q[k + 1][j]) < accuracy)
                Q[k + 1][j] = 0.0;
        }
        for (std::size_t j = k; j < cols; j++){
            Type temp = matrix[k][j];
            matrix[k][j] = c * matrix[k][j] + s * matrix[k + 1][j];
            matrix[k + 1][j] = -s * temp + c * matrix[k + 1][j];
            if (std::abs(matrix[k + 1][j]) < accuracy)
                matrix[k + 1][j] = 0.0;
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
std::size_t transposeMatrix(std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return 0;
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = i + 1; j < cols; j++){
            Type temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
    return rows;
}

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    std::vector<Type> b1(rows);  // Правая часть после подстановки полученного решения
        for (std::size_t i = 0; i < rows; i++){
                Type sum = 0.0;
                for (std::size_t k = 0; k < cols; k++){
                    sum += lCoefs[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(rows); // Вектор невязки
        for (std::size_t i = 0; i < rows; i++){
            discrepancyVector[i] = rCoefs[i] - b1[i];
        }
        Type discrepancy = 0.0; // Невязка
        for (std::size_t i = 0; i < rows; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
    return std::sqrt(discrepancy);
}

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type norm1OfMatrix = 0;
    for (std::size_t j = 0; j < cols; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < rows; i++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > norm1OfMatrix)
            norm1OfMatrix = sum;
    }
    return norm1OfMatrix;
}

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type norm1OfMatrix = findMatrixNorm1(matrix);
    Type norm1OfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        norm1OfInvMatrix = findMatrixNorm1(invMatrix);
    else
        norm1OfInvMatrix = INFINITY;
    Type cond = norm1OfMatrix * norm1OfInvMatrix;
    return cond;
}

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type normInfOfMatrix = 0;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > normInfOfMatrix)
            normInfOfMatrix = sum;
    }
    return normInfOfMatrix;
}

template<typename Type>
Type normOfMatrix(const std::vector<std::vector<Type>> &matrix, double p){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (p == 1.0){
        return findMatrixNorm1(matrix);
    }
    if (p == INFINITY){
        return findMatrixNormInf(matrix);
    }
    return NAN;
}

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type normInfOfMatrix = findMatrixNormInf(matrix);
    Type normInfOfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        normInfOfInvMatrix = findMatrixNormInf(invMatrix);
    else
        normInfOfInvMatrix = INFINITY;
    Type cond = normInfOfMatrix * normInfOfInvMatrix;
    return cond;
}

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix,
    SOLUTION_FLAG (*method)(std::vector<std::vector<Type>> &, std::vector<Type>&, std::vector<Type>&, Type accuracy)){
    std::size_t rows = inputMatrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = inputMatrix[0].size();
    else
        return NOT_INVERTIBLE;
    if (rows != cols)
        return NOT_INVERTIBLE;
    resMatrix.resize(rows);
    std::vector<std::vector<Type>> E(rows);
    for (std::size_t i = 0; i < rows; i++){
        E[i].resize(cols, 0);
        resMatrix[i].resize(cols);
    }
    for (std::size_t i = 0; i < rows; i++){
        E[i][i] = 1;
    }
    std::vector<Type> solution(rows);
    std::vector<std::vector<Type>> tempMatrix(rows);
    SOLUTION_FLAG flag;
    for (std::size_t i = 0; i < rows; i++){
        tempMatrix = inputMatrix;
        flag = (*method)(tempMatrix, E[i], solution, 1e-7);
        if (flag == NO_SOLUTION){
            for (std::size_t i = 0; i < rows; i++)
                resMatrix[i].clear();
            resMatrix.clear();
            return NOT_INVERTIBLE;
        }
        for (std::size_t k = 0; k < rows; k++)
            resMatrix[k][i] = solution[k];
    }
    return IS_INVERTIBLE;
}

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++)
        sum += std::abs(vector[i]);
    return sum;
}

template<typename Type>
Type norm2OfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++){
        sum += std::pow(vector[i], 2);
    }
    return std::sqrt(sum);
}

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type max = std::abs(vector[0]);
    for (std::size_t i = 1; i < vector.size(); i++)
        if (std::abs(vector[i]) > max)
            max = std::abs(vector[i]);
    return max;
}

template<typename Type>
Type normOfVector(const std::vector<Type> &vector, double p){
    if (!vector.size()){
        return NAN;
    }
    if (p == 2.0){
        return norm2OfVector(vector);
    }
    if (p == 1.0){
        return norm1OfVector(vector);
    }
    if (p == INFINITY){
        return normInfOfVector(vector);
    }
    return NAN;
}

template<typename Type>
Type findLowerBoundOfcond1(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
Type findLowerBoundOfcondInf(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec, std::vector<Type> &result){
    std::size_t rows1 = matrix.size();
    std::size_t cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return NOT_MULTIPLIED;
    std::size_t rows2 = vec.size();
    if (cols != rows2)
        return NOT_MULTIPLIED;
    result.resize(rows1);
    std::size_t rows = result.size();
    for (size_t i = 0; i < rows; i++){
        Type sum = 0;
        for (size_t k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return IS_MULTIPLIED;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else{
        os << 0;
        return os;
    }
    for (std::size_t i = 0; i < rows - 1; i++){
        for (std::size_t j = 0; j < cols; j++){
            os << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
    for (std::size_t j = 0; j < cols; j++){
            os << matrix[rows - 1][j] << ' ';
        }
    return os;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<Type> &vector){
    std::size_t rows = vector.size();
    if (!rows){
        os << 0;
        return os;
    }  
    os << "{ ";
    for (std::size_t i = 0; i < rows - 1; i++)
        os << vector[i] << ", ";
    os << vector[rows - 1] << ' ';
    os << '}';
    return os;
}

template<typename Type>
std::vector<Type> operator+(const std::vector<Type>& vec1, const std::vector<Type>& vec2)
{
    auto result = vec1;
    const std::size_t size_max = std::max<std::size_t>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (std::size_t i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] += vec2[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator-(const std::vector<Type>& vec1, const std::vector<Type>& vec2)
{
    auto result = vec1;
    const std::size_t size_max = std::max<std::size_t>(vec1.size(), vec2.size());
    result.resize(size_max, 0);
    for (std::size_t i = 0; i < result.size() && i < vec2.size(); i++){   
        result[i] -= vec2[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec)
{
    std::size_t rows1 = matrix.size();
    std::size_t cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return std::vector<double>(1, NAN);
    std::size_t rows2 = vec.size();
    if (cols != rows2)
        return std::vector<double>(rows2, NAN);
    std::vector<Type> result(rows1);
    for (std::size_t i = 0; i < rows1; i++){
        Type sum = 0;
        for (std::size_t k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(Type num, const std::vector<Type> &vec){
    std::size_t size = vec.size();
    if (!size)
        return std::vector<double>(1, NAN);
    std::vector<Type> res(size);
    for (std::size_t i = 0; i < size; i++){
        res[i] = num * vec[i];
    }
    return res;
}

template<typename Type>
std::vector<std::vector<Type>> operator*(Type num, const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (!rows)
        return std::vector<std::vector<Type>>(1, std::vector<Type>(1, NAN));
    else
        cols = matrix[0].size();
    std::vector<std::vector<Type>> res;
    std::vector<Type> tempVec(cols);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++)
        {
            tempVec[j] = num * matrix[i][j];
        }
        res.push_back(tempVec);
    }
    return res;
}

template<typename Type>
std::size_t simpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type tao, Type accuracy, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t JacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t relaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, Type omega, double p, Type epsilon_0, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t relaxationMethodFor3Diag(const std::vector<Type> &a, const std::vector<Type> &b, const std::vector<Type> &c, const std::vector<Type> &d, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, Type accuracy, Type omega, double p, Type epsilon_0, std::size_t stopIt){
    if (!b.size() || b.size() != d.size() || a.size() != b.size() - 1 || c.size() != a.size())
        return NO_SOLUTION;
    std::size_t dim = b.size();
    solution.resize(dim); // Искомое решение
    std::vector<Type> prev_solution = firstVec;
    solution[0] = (1 - omega) * prev_solution[0] - (omega / b[0]) * (c[0] * prev_solution[1] - d[0]);
    for (std::size_t i = 1; i < dim - 1; i++){    
        solution[i] = (1 - omega) * prev_solution[i] - (omega / b[i]) * (a[i - 1] * solution[i - 1] + c[i] * prev_solution[i + 1] - d[i]);
    }
    solution[dim - 1] = (1 - omega) * prev_solution[dim - 1] - (omega / b[dim - 1]) * (a[dim - 2] * prev_solution[dim - 2] - d[dim - 1]);
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (diffNorm / (normOfVector(prev_solution, p) + epsilon_0) > accuracy){
        prev_solution = solution;
        solution[0] = (1 - omega) * prev_solution[0] - (omega / b[0]) * (c[0] * prev_solution[1] - d[0]);
        for (std::size_t i = 1; i < dim - 1; i++){    
            solution[i] = (1 - omega) * prev_solution[i] - (omega / b[i]) * (a[i - 1] * solution[i - 1] + c[i] * prev_solution[i + 1] - d[i]);
        }
        solution[dim - 1] = (1 - omega) * prev_solution[dim - 1] - (omega / b[dim - 1]) * (a[dim - 2] * prev_solution[dim - 2] - d[dim - 1]);
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

// x = C * x + y 
template<typename Type>
QUADRATIC_FLAG findCanonicalFormSimpleIt(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type tao){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    } 
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                C[i][j] = -tao * lCoefs[i][j];
            }
            else{
                C[i][j] = -tao * lCoefs[i][j] + 1;
            }
        }
    }
    for (std::size_t i = 0; i < cols; i++){
        y[i] = tao * rCoefs[i];
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormJacobi(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    } 
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                C[i][j] = -lCoefs[i][j] / lCoefs[i][i];
            }
            else{
                C[i][j] = 0.0;
            }
        }
    }
    for (std::size_t i = 0; i < rows; i++){
        y[i] = rCoefs[i] / lCoefs[i][i];
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    }
    std::vector<Type> columnOfMatrix(rows, 0);
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * C[j][k];
            }
            Type sum2 = 0.0;
            if (k > i){
                sum2 = lCoefs[i][k];  
                C[i][k] = -(omega / lCoefs[i][i]) * (sum1 + sum2);
            }
            else{
                if (k < i){
                    C[i][k] = -(omega / lCoefs[i][i]) * sum1;
                }
                else{
                    C[i][k] = (1 - omega) - (omega / lCoefs[i][i]) * sum1;
                }
            }
        }
    }
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        Type sum = omega * rCoefs[i] / lCoefs[i][i];
        for (std::size_t k = 0; k < i; k++){
            sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * y[k]; 
        }
        y[i] = sum;
    }
    return IS_QUADRATIC;
}

template<typename Type>
QUADRATIC_FLAG findCanonicalFormRelaxation2(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, std::vector<std::vector<Type>> &C, std::vector<Type> &y, Type omega){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    cols = rCoefs.size();
    if (rows != cols)
        return NOT_QUADRATIC;
    C.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        C[i].resize(cols, 0);
    }
    C[0][0] = 1 - omega; 
    for (std::size_t j = 1; j < cols; j++){
        C[0][j] = -omega * lCoefs[0][j] / lCoefs[0][0];
    }
    for (std::size_t i = 1; i < rows; i++){
        for (std::size_t j = 0; j < i; j++){
            Type sum = 0.0;
            for (std::size_t k = 0; k < i; k++){
                sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * C[k][j];
            }
            C[i][j] = sum;
        }
        for (std::size_t j = i; j < cols; j++){
            Type sum = 0.0;
            if (i != j){
                sum = -omega * lCoefs[i][j] / lCoefs[i][i];
            }
            else{
                sum = 1 - omega;
            }
            for (std::size_t k = 0; k < i; k++){
                sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * C[k][j];
            }
            C[i][j] = sum;
        }
    }
    y.resize(cols, 0);
    for (std::size_t i = 0; i < rows; i++){
        Type sum = omega * rCoefs[i] / lCoefs[i][i];
        for (std::size_t k = 0; k < i; k++){
            sum += -omega * (lCoefs[i][k] / lCoefs[i][i]) * y[k]; 
        }
        y[i] = sum;
    }
    return IS_QUADRATIC;
}

template<typename Type>
Type findLowerBoundOfIterations(const std::vector<std::vector<Type>> &lCoefs, 
const std::vector<Type> &rCoefs, const std::vector<Type> &firstVec, Type accuracy, ITERATION_METHOD_FLAG method, Type tao, Type omega, double p){
    std::size_t rows = lCoefs.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0.0;
    std::vector<Type> prev_solution = firstVec;
    std::vector<Type> solution(rows);
    std::vector<std::vector<Type>> C;
    std::vector<Type> y;
    switch (method){
        case SIMPLE_IT:
            for (std::size_t i = 0; i < rows; i++){
                Type sum = 0.0;
                for (std::size_t j = 0; j < cols; j++){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
                solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
            }
            findCanonicalFormSimpleIt(lCoefs, rCoefs, C, y, tao);
            break;

        case JACOBI:
            for (std::size_t i = 0; i < rows; i++){
                    Type sum = 0.0;
                    for (std::size_t j = 0; j < cols; j++){
                        if (i != j){
                        sum += lCoefs[i][j] * prev_solution[j];
                        }
                    }
                solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
            }
            findCanonicalFormJacobi(lCoefs, rCoefs, C, y);
            break;

        case RELAXATION:
            for (std::size_t i = 0; i < rows; i++){
                Type sum1 = 0.0;
                for (std::size_t j = 0; j < i; j++){
                    sum1 += lCoefs[i][j] * solution[j];
                }
                Type sum2 = 0.0;
                for (std::size_t j = i + 1; j < cols; j++){
                    sum2 += lCoefs[i][j] * solution[j];
                }
                solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
            }
            findCanonicalFormRelaxation(lCoefs, rCoefs, C, y, omega); 
            break;
            default:
                break;
    }
    Type distance = normOfVector(solution - prev_solution, p);
    Type q = normOfMatrix(C, p);
    return std::ceil(std::log(accuracy * (1 - q) / distance) / std::log(q));
}

template<typename Type>
std::size_t findExactItersSimpleItMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type tao, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t findExactItersJacobiMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
std::size_t findExactRelaxationMethod(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type accuracy, Type omega, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (normOfVector(solution - rightSolution, p) > accuracy){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return numOfIt;
}

template<typename Type>
Type findNormOfErrAfterEstIt_SIT(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type tao, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += lCoefs[i][j] * prev_solution[j];
        }
        solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                sum += lCoefs[i][j] * prev_solution[j];
            }
            solution[i] = prev_solution[i] + tao * (rCoefs[i] - sum);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

template<typename Type>
Type findNormOfErrAfterEstIt_JAC(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return INFINITY;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            if (i != j){
                sum += lCoefs[i][j] * prev_solution[j];
            }
        }
        solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum = 0.0;
            for (std::size_t j = 0; j < cols; j++){
                if (i != j){
                    sum += lCoefs[i][j] * prev_solution[j];
                }
            }
            solution[i] = (1/lCoefs[i][i]) * (rCoefs[i] - sum);
        }   
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt)
            break;
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

template<typename Type>
Type findNormOfErrAfterEstIt_REL(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, 
const std::vector<Type> &firstVec, std::vector<Type> &solution, std::vector<Type> &rightSolution, Type bound, Type omega, Type accuracy, double p, std::size_t stopIt){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return 0;
    std::vector<Type> prev_solution = firstVec;
    for (std::size_t i = 0; i < rows; i++){
        Type sum1 = 0.0;
        for (std::size_t j = 0; j < i; j++){
            sum1 += lCoefs[i][j] * solution[j];
        }
        Type sum2 = 0.0;
        for (std::size_t j = i + 1; j < cols; j++){
            sum2 += lCoefs[i][j] * solution[j];
        }
        solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
    }
    std::vector<Type> diffVec = solution - prev_solution;
    Type diffNorm = normOfVector(diffVec, p);
    std::size_t numOfIt = 1;
    while (numOfIt < bound + 1){
        prev_solution = solution;
        for (std::size_t i = 0; i < rows; i++){
            Type sum1 = 0.0;
            for (std::size_t j = 0; j < i; j++){
                sum1 += lCoefs[i][j] * solution[j];
            }
            Type sum2 = 0.0;
            for (std::size_t j = i + 1; j < cols; j++){
                sum2 += lCoefs[i][j] * solution[j];
            }
            solution[i] = (1 - omega) * prev_solution[i] - (omega / lCoefs[i][i]) * (sum1 + sum2 - rCoefs[i]);
        }
        diffVec = solution - prev_solution;
        diffNorm = normOfVector(diffVec, p);
        if (numOfIt == stopIt){
            break;
        }
        numOfIt++;
    }
    return normOfVector(solution - rightSolution, p);
}

// Лаба 3
template<typename Type>
std::size_t findEigenNumsQRMethod(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy, bool hasShift){
    std::size_t numOfIters = 0; // Количество итераций
    Type shift = 0.0; // Сдвиг
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return numOfIters; //////////////////
    eigenList.resize(rows, 0.0);
    std::size_t eigenRow = rows - 1; // Строка в которой по итогу вращений должно получиться собственное значение
    std::vector<std::vector<Type>> Q;
    std::vector<std::vector<Type>> R;
    R.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        R[i].resize(cols, 0.0);
    }
    ///////////////////////////////////////////////////////////////////////////////
    std::string IN_FILE_PATH_4 = "D:\\Calc_Methods\\Lab3\\Tests\\test4.txt";
    std::string IN_FILE_PATH_5 = "D:\\Calc_Methods\\Lab3\\Tests\\test5.txt";
    std::string IN_FILE_PATH_6 = "D:\\Calc_Methods\\Lab3\\Tests\\test6.txt";
    ///////////////////////////////////////////////////////////////////////////////
    while (eigenRow != 0){
        numOfIters++;
        ///////////////////////////////////////////////////        
        writeMatrixFile(matrix, IN_FILE_PATH_4);
        ///////////////////////////////////////////////////
        // Сдвиг
        if (hasShift){
            shift = matrix[eigenRow][eigenRow];
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] -= shift; 
        }
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                R[i][j] = matrix[i][j];
            }
        }
        findQMatrix(R, Q, accuracy);
        ///////////////////////////////////////////////////////////
        writeMatrixFile(Q, IN_FILE_PATH_5);
        writeMatrixFile(R, IN_FILE_PATH_6);
        //////////////////////////////////////////////////////////
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                Type sum = 0.0;
                for (std::size_t k = i; k < rows; k++){
                    sum += R[i][k] * Q[k][j];
                }
                matrix[i][j] = sum;
            }
        }
        // Обратный сдвиг
        if (hasShift){
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] += shift; 
            }
        }
        Type sumOfEigenRow = 0.0;
        for (std::size_t j = 0; j < eigenRow; j++){
            sumOfEigenRow += std::abs(matrix[eigenRow][j]);
        }
        if (sumOfEigenRow < accuracy){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow]; 
            for (std::size_t i = 0; i < eigenRow; i++){
                Q[i].resize(eigenRow);
                R[i].resize(eigenRow);
            }
            Q.resize(eigenRow);
            R.resize(eigenRow);
            eigenRow--;
            rows--;
            cols--;
        }
        if (eigenRow == 0){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow];
            break;
        }
    }
    return numOfIters;
}

template<typename Type>
QUADRATIC_FLAG getHessenbergMatrix(std::vector<std::vector<Type>> &matrix, Type accuracy, bool isSymmetric){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NOT_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    if (!isSymmetric){
        for (std::size_t k = 1; k < rows; k++){
            for (std::size_t i = k + 1; i < rows; i++){
                if (std::abs(matrix[i][k - 1]) >= accuracy){
                    Type c = matrix[k][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    Type s = matrix[i][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    for (std::size_t j = k - 1; j < cols; j++){
                        Type temp = matrix[k][j];
                        matrix[k][j] = c * matrix[k][j] + s * matrix[i][j];
                        matrix[i][j] = -s * temp + c * matrix[i][j];
                    }
                    for (std::size_t j = k - 1; j < rows; j++){
                        Type temp = matrix[j][k];
                        matrix[j][k] = c * matrix[j][k] + s * matrix[j][i];
                        matrix[j][i] = -s * temp + c * matrix[j][i];
                    }
                }
            }
        }    
    }
    else{
       for (std::size_t k = 1; k < rows; k++){
            for (std::size_t i = k + 1; i < rows; i++){
                if (std::abs(matrix[i][k - 1]) >= accuracy){
                    Type c = matrix[k][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    Type s = matrix[i][k - 1] / std::sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[i][k - 1] * matrix[i][k - 1]);
                    for (std::size_t j = k - 1; j < cols; j++){
                        Type temp = matrix[k][j];
                        matrix[k][j] = c * matrix[k][j] + s * matrix[i][j];
                        matrix[j][k] = matrix[k][j];
                        matrix[i][j] = -s * temp + c * matrix[i][j];
                        matrix[j][i] = matrix[i][j];
                    }
                }
            }
        } 
    }
    return IS_QUADRATIC;
}

template<typename Type>
std::size_t findEigenNumsQRMethodHessenberg(std::vector<std::vector<Type>> &matrix, std::vector<Type> &eigenList, Type accuracy, bool hasShift){
    std::size_t numOfIters = 0; // Количество итераций
    Type shift = 0.0; // Сдвиг
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return numOfIters; 
    eigenList.resize(rows, 0.0);
    std::size_t eigenRow = rows - 1; // Строка в которой по итогу вращений должно получиться собственное значение
    std::vector<std::vector<Type>> Q;
    std::vector<std::vector<Type>> R;
    R.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        R[i].resize(cols, 0.0);
    }
    // Приводим к матрице Хессенберга
    getHessenbergMatrix(matrix, accuracy);
    while (eigenRow != 0){
        numOfIters++;
        // Сдвиг
        if (hasShift){
            shift = matrix[eigenRow][eigenRow];
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] -= shift; 
        }
        }
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                R[i][j] = matrix[i][j];
            }
        }
        findQMatrix3Diag(R, Q, accuracy);
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = 0; j < cols; j++){
                Type sum = 0.0;
                for (std::size_t k = i; k < rows; k++){
                    sum += R[i][k] * Q[k][j];
                }
                matrix[i][j] = sum;
            }
        }
        // Обратный сдвиг
        if (hasShift){
            for (std::size_t i = 0; i < rows; i++){
                matrix[i][i] += shift; 
            }
        }
        Type sumOfEigenRow = 0.0;
        for (std::size_t j = 0; j < eigenRow; j++){
            sumOfEigenRow += std::abs(matrix[eigenRow][j]);
        }
        if (sumOfEigenRow < accuracy){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow]; 
            for (std::size_t i = 0; i < eigenRow; i++){
                Q[i].resize(eigenRow);
                R[i].resize(eigenRow);
            }
            Q.resize(eigenRow);
            R.resize(eigenRow);
            eigenRow--;
            rows--;
            cols--;
        }
        if (eigenRow == 0){
            eigenList[eigenRow] = matrix[eigenRow][eigenRow];
            break;
        }
    }
    return numOfIters;
}