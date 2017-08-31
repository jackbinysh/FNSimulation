function [ totalerror ] = CurveOffsetError(offset, Curvatures1,IntegratedLengths1,Curvatures2,IntegratedLengths2 )
%CURVEOFFSETERROR the error between two curves offset from one another

% lets construct an extended version of our periodic curve
SecondCurveLength = IntegratedLengths2(end);
ExtendedLengths2 = vertcat(IntegratedLengths2-SecondCurveLength, IntegratedLengths2, IntegratedLengths2 + SecondCurveLength);
ExtendedCurvatures2 = vertcat(Curvatures2, Curvatures2, Curvatures2);

%errors = Curvatures1 - interp1(ExtendedLengths2,ExtendedCurvatures2, IntegratedLengths1 - offset);
%totalerror = sum(abs(errors));

totalerror = arrayfun( @(i)(sum(abs(Curvatures1 - interp1(ExtendedLengths2,ExtendedCurvatures2, IntegratedLengths1 - offset(i))))), 1:length(offset) );



end

