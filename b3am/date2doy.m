function doy = date2doy(indate)

d1 = datevec(indate);
d1 = d1(1:3);
d2 = d1;
d2(2:3) = 0;

doy = datenum(d1)-datenum(d2);

end