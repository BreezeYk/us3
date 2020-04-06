function output= mutualExchange(signature)
    p2=size(signature,2);
    p3=size(signature,1);
    newSignature=zeros(p3,p2);
    newSignature (:,:,1:p3/2) = signature(:,:,p3/2+1:p3);
    newSignature (:,:,p3/2+1:p3) = signature(:,:,1:p3/2);
    output= newSignature;