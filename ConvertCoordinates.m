function [Xh,Yh] = ConvertCoordinates(axis, Xa,Ya)
XLim     = axis.XLim; YLim = axis.YLim; 
Position = axis.Position;
Xh       = Position(1) + Position(3)*(Xa - XLim(1))/diff(XLim);
Yh       = Position(2) + Position(4)*(Ya - YLim(1))/diff(YLim);