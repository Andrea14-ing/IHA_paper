function H0 = Hinv(H)
%inverse of a RotoTranslation matrix
H0=H;
H0(1:3,1:3,:)=permute(H0(1:3,1:3,:),[2,1,3]);
o=-sum(H(1:3,1:3,:).*H(1:3,[4 4 4],:));
H0(1:3,4,:)=permute(o,[2,1,3]);
end
