%% MFN 1 (Chapter 10)

% 10.2; Vec. Addition and Multiplication
% 3 types of operations applicable to vectors
    % Vec. addition
        % summing two vecs. element by element (dims must match)
    % Scalar multiplication
        % scalar= single #; so scalar multiplication means multiplying each
        % element in the vector by a single number
            % this only renders stretching or shrinking a vec rather than a
            % directional change 
    % Vec. multiplication
        % has 2 forms,
            % inner product
                % AKA the dot. product, the dot product is a sigle number
                % that tells you something about the relationship between two vectors,
               
                % row vec (1x3) * column vec (3x1) = dot prod. (1x1)
                % i.e. 
                    % v = [3 2 6];  w = [-1 3 2];
                    % v.w = [3*-1 + 2*3 + 6*2] = 15

            % outer product
                % is a rank 1 matrix, not important, but important to know
                % that for this book, if you get an outer product, its an
                % error

                % column vec (3x1) * row vec (1x3) = outer prod. (3x3)
    % Vector Length 
        % the length of a vector can be thought of as the hypotenuse of a right
        % triangle. in this case the length is the sqr root of the sum of the
        % two sides squared. for vectors, we can use the sqr root of the dot
        % product of the *vector with itself*
    % The Dot Product Geometrically 
        % geometrically the dot product is the multiplication of the lgnths
        % of two vectors, scaled by the cosine of the angle between them 
            % the formula for geometric interpretation 
                % = |v||w|cos θ
            % the smaller the angle between two vectors, the larger the dot
            % product
    % Ways to Compute the Dot Prod 
        % sum(v.*w)
        % dot(v,w)
        % v'*w

% 10.3 Matrices
    % matrices can be though of as a stack of row vectors, or a row of
    % column vectors
        % the form is always M x N, Rws x cols
    % Special Matrices 
        % Square- shaped like a square, or, have same number of rows and
        % cols 
            % A non-square matrix is a rectangular matrix
    % Symmetric Matrix 
        % a square matrix where the upper triangle is the same as the
        % lower, essentially a mirror image across the diagonal, and if
        % transposed ('), the matrix equals itself, in value not size
    % Identity Matrix
        % the main diagonal is all 1's, all other values are zeros
            % the identity matrix, times any other matrix or vector, ='s
            % that same vector. Its like multiplying by 1
                % syntax for this is 'eye(N)'
% 10.4 Finding Your Way around a Matrix    
    % matrix locations are encoded as row,column. V(3,5)
        % multiple elements: V(3:8,10:20)
        % noncontigious elements: V([1 3 9],[10 20 30])
        % WARNING: when indexing always make sure there are V-1 commas in your code
   % There are two matrix formats/organizations
        % Row 1 at the bottom: 'xy-format'
            % used mostly for images
        % Row 1 at the top: 'ij-format'
            % used mostly for matrices 
% 10.5 Matrix Multiplication 
    % not all matrices may be multiplied, ONLY those with matching inner
    % dimensions!
        % the size outcome of this multiplication is the outer dimensions
    % Matrix multiplication is Non-Commutative, 
        % meaning AB is not the same
        % as BA, this makes sense, as the inner dimensions much match, and the
        % product = the outer dimensions
    % Matrices that are multipliable are called conformable 
    % Matrix multiplication is simply the dot prod of each ROW of the left
    % matrix, with each COLUMN of the right matrix
    % Transposesand matrix multiplication
        % the transpose of a matrix multiplication is the multiplication of 
        % the individual matrices transposed and in reverse order.
            % translated into matlab code 
                % (A*B)' = (B')*(A')
% 10.6 When to use .* or ./ vs * and /
    % the difference comes from the differences in multiplication for
    % scalars versus vectors and matrices
        % the * is for matrix multiplication 
        % * or .* can be used for scalars (with anything else)
    % However, .* is for scalar multiplication, ie if used on a matrix .*
    % would yield element-wise multiplication
        % eg if A and B are both 3x4 matrices, A*B would be an err (because
        % matrices multiplied together must have homogenous inside values
        % (3x4 3x4)
    % Confusion can arise with equal sized matrices
        % eg A=4x4, B=4x4, both A*B and A.*B = 4x4, however the two are
        % vastly different, as standard practice .* should be the default
        % syntax for multiplication, and * only used when youre sure of
        % matrix multiplication intent (rmember matrix multiplication = the
        % dot prod of every row and every column, and equals in size, the
        % inner dimensions of the two things being multiplied)

% 10.7 Linear Independence and Rank 
    % linear independence is a property of a set of vectors (individual or
    % in the context of columns/rows of a matrix)
    
    % The algebraic interpretation is that if one vector can be produced
    % from a linear combination (scaling and adding of other vectors), then that
    % set of vectors are considered linearly depedent 
        %eg a=[4 5], e=[2 8], and f=[-11 -11]...
        % linearly combining vectors a and e (-3a + .5e) will equal f, and
        % thus f provides no novel information
    
    % Geometrically, a set of vectors is linearly independent if each
    % vector points in a diff geometric dimension. 
        % eg vectors [1 2] and [2 4] are not independent, they point in the
        % same direction 
            % but vectors [1 2], and [2 5] point in diff directions, and
            % thus would be dependent 
                % a third vector would be dependent by virtue of the fact
                % that the first 2 vectors cover the 2 dimensions of a 2D
                % space
                    % Further; any third vector in 2Dimensions can always
                    % be created by adding and subtracting scaled versions
                    % of the other two independent vectors 
    % The strongest case of linear independence is orthogonality.
    % Algebraically, two vectors are orthogonal when the dot product
    % between them is zero. Geometrically two vectors meet at a 90 deg.
    % angle 
    
    % Rank refers to the number of linearly independent columns/rows in a matrix
        % if every row/col in a matrix formed a linearly independent set of
        % vecs, that matrix is considered 'full rank'
            % matlab syntax for rank calculation is 'rank(A)'

% 10.8 The Matrix Inverse
    % to neutralize a scalar you can multiply it by its inverse,
    % yielding a 1 (e.g. 1/5 * 5)

    % Analogous to this (multiplier), is the matrix inverse, written A^-1
    % However, we must remember that the neutral matrix (ident matrix), 
    % is not all 1's
    % but as previously introduced, (ones on diagonal, zeros everywhere
    % else)
        % so A * A^-1 = Ident. matrix A
    
    % There are quite a few stipulations on matrix inversion 
        % only square matrices can be inverted 
            % within those, only one's with full rank may be inverted
    % non invertible matrices are called 'singular' 

    % a pseudo inverse can be used, a close approximation for a true
    % inverse and equals the true inverse when the matrix is squared
    % and invertible (Q: i thought the whole point is its not
    % invertible in the first place here?)
 
% 10.9 Solving Ax = b
    % A MAIN function of using matrix inverse is moving a matrix from
    % one side of an equation to the other, ie in the case of solving
    % systems of linear equations and, least-squares solutions
        % eg (with scalars), solve for x. ax = b
            % obvi x = b/a is the solution
        % its the same for matrices; multiplying both sides by a^-1
    % least-squares equations for solving systems of linear equations
        % the least squares approach is really important
        % the basic idea is solve Ax = b
            % A is a matrix of predictors or indepenedent variables
                % (columns being variables, rows as trials)
            % x is a vector containing the coeffs (or regression weights)
            % b is a vector of observed data points 
        % So, Ax = b is asking the question: what linear combination of
        % independent variables can explain the observed data?
            % the goal here is to find the vector of unknowns, or x...
                % and to solve we simply must move A to the other sider of
                % the equation, using the matrix inverse 
                    % playsout like this. 
                        % Ax = b
                        % (A^-1A)x = A^-1b
                        % Ix = A^-1b
                        % and A^-1 is NxN, while b is Nx1 so it works
                        % out
%10.10 Making Symmetric Squares from Rectangles
    % Typically however, A will not be invertible, because usually it wont
    % be a square (specially given its form of variables x trials)
        % So we need to make it a square matrix so that it can be inverted
        % and solved for???
    % the solution is to use the transpose, 
        % luckily a matrix times its
        % transpose is always a square matrix, obvi cause NxM matrix * MxN
        % matrix = NxN matrix 
    % additionally a matrix times its transpose, is guaranteed to be
    % symmetric 
        % square symmetric matrices have great properties, sucha as
        % orthogonal eigenvectors 
    % For the least squares equation the most important thing is that
    % square symmetric matrices are invertible as long as the columns of A
    % are linearly independent (AKA the rank of the matrix is equal to the
    % number of columns)

    % so to solve the least-squares equation (Ax = b) you start by
    % multiplying both sides by A^T then proceeding.
    % as in the previous section
        % Ax = b 
        % A^T(A)x = b(A^T)
        % A^TAx = A^Tb
        % (A^TA)^-1 (A^TA)x= (A^TA)^-1 (A^T)b
        % x = (A^T A)^-1 (A^T) b
    % this solution has a dedicated representation in MATLAB, as x =
    % (A'*A)\A'*b;
%10.11 Full and Sparse Matrices
    % a sparse matrix is one where a majority of the elements are zeros and
    % is represented in matlab in a unique way, by the coordinates of the
    % non-zero value and the value 
    

% Excersises

% 1a. valid, 4x4, 1b.invalid, 1c. invalid, d. valid, 4x5, e. invalid, f.
% valid, 4x4, g. valid, 5x4, h. invalid
    % Test 
        A = rand(4,4); B = rand(4,7); C = rand(5,7);
        A*A; B*A; C*A; B*C'; C*B*A; B*B'; C*B'*A; C*C;
% 2 
    a = [1,-3]; b = [3,1]
    dot(a,b)
    a = [1,-5]; b = [3,1]
    % no, -3 is the only solution for orthogonality between these vectors
    a = [1,-5]; b = [5,1]
% 3 
    aNumber = randn; round(1000*aNumber)/1000;
% 4
    % reset
        format;
% 5
    %a
    plot3([0 v3(1)],[0 v3(2)])
    plot3([0 v3(1)],[0 v3(2)],[0 v3(3)])    
    v3= [6 -4 5 1];figure; plot([0,v3(1)]); hold on; plot([0, v3(2)])
    %b
    a = rand(2,3,4); a'
    a = rand(2,3,4); a(:,:,1)'
    %c
    inv(rand(2,3))
    inv(rand(3,3))
 % 6
    a = [2,4;4,8];
    ainverse = inv(a);
    ainverse = pinv(a);
    a3 = a*ans

    a = [1,5;1,-5];
    rank(a)
    ainverse = inv(a);
    ainverse = pinv(a);
    a3 = a*ainverse   
    % No, for a full rank matrix, if you multiply by the pinverse you get
    % something else, if you use the normal inverse (which works without
    % error), you get the identity matrix

 % 7
    a = rand(2,3)
    % AA^T = 2x2
    % A^T (A) = 3x3
    first = (a)*(a')
    plot([0 first(1,1)],[0 first(1,2)],[0 first(2,1)], [0 first(2,2)]);
    % that they're on the same line makes them symmetric
    second = (a')*(a)
    plot([0 first(1,1)],[0 first(1,2)],[0 first(2,1)], [0 first(2,2)]);
    axis([0 2 0 2])

 % 8 
     (4,2) or (2,4) %how can you tell past that?? I mean i know a symmetric
     %matrix should be reflective across the diagonal, but either of those
     %two numbers could be wrong??

 % 9
    a = rand(1,3,2);
    SuDot = sum(a(:,:,1).*a(:,:,2))
    realDot = dot(a(:,:,1),a(:,:,2))
    aneg = -1.*(a)
    SuDot = sum(aneg(:,:,1).*aneg(:,:,2))
    % the dot product can be negative if you're computing it between a
    % positive and negative value...I think, ofc w/the dot of a vector and
    % itself, itll never be negative cause (-)(-) = + 

  % 10 
     a = [1 3 3; 0 -1 -2; 7 2 -17];
     rank = rank(a);
     % i dont think I'm suprised, i guess cause the column/row with the '0'
     % and either of the other rows/columns will be independent. Thus, I
     % would change the zero to change rank
     a2 = [1 3 3; 4 -1 -2; 7 2 -17];
     rank = rank(a2);

  % 11
     a = rand(4,8)
     x = 1; %yourRow
     y = 3; %yourCol
     lindex = (x)+(length(a(:,1))*(y-1)) 
     % to get the linear index you take the row input + # rows in the matrix * how far you're from the 1st row
     % check
     sub2ind([4,8],1,3);

   % 12
     a = rand(1,4)
     b = rand(3,4)
     c = repmat(a,3,1)
     %rank should be zero; wrong, 1; oh misconception, I guess 1 refers to
     %1 independent vector being present
     rank(b)
     pointwise = b.*c

   % 13
     bettersol = bsxfun(@times,b,a) % this works, tbh i was a bit confused 
     % about what function handle to use, we wanted point-wise
     % multiplication, not array multiplication?? (I guess array insinuates
     % point wise; Cause it's not matrix multiplication...)
   
   % 14
     diag([2 1 5 6])
     a = diag([1 1 1 1 1 1 1 1 1 1]) % using only diag

     b= [1]
     b = repmat(b,1,10)
     b = diag(b) % using diag and repmat

   % 15
     a = [1 3 2]
     b = [5 2 3]
     dot(a,b) % i was right in my hand computation
     % trying to remember the mat-prod. notation w.out notes (x^Ty);
     % it would be x':3x1 *(matrixmult) y:1x3 (this looks good! == 1x1 aka
     % dot product)
     dotBynotation = a*b' % coming out 3x3... seems like its doing point-wise mult instead of matrix
     % Nope! silly mistake, 1x3 x 1x3, i need to transpose the 2nd vector
     % to get 1x1, the above is now correct.




      
