function y = quat2euler(x)

    q = [x(1) -x(2) -x(3) -x(4)];
    R(1,1,:) = 2.*q(:,1).^2-1+2.*q(:,2).^2;
    R(2,1,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
    R(3,1,:) = 2.*(q(:,2).*q(:,4)+q(:,1).*q(:,3));
    R(3,2,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
    R(3,3,:) = 2.*q(:,1).^2-1+2.*q(:,4).^2;

    phi = atan2(R(3,2,:), R(3,3,:) );
    theta = -atan(R(3,1,:) ./ sqrt(1-R(3,1,:).^2) );
    psi = atan2(R(2,1,:), R(1,1,:) );

    y = [real(phi(1,:)') real(theta(1,:)') real(psi(1,:)')]';
    
end

