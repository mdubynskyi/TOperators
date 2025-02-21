
#include <vector>

class SO6{
    public:
        SO6();
        // SO6(std::vector<int8_t> t);
        // SO6(Z2[6][6], std::vector<int8_t> t); //initializes matrix according to a 6x6 array of Z2
        SO6 operator*(SO6&); //mutliplication
        void fixSign();
        void lexOrder();
        inline Z2& operator()(int8_t col, int8_t row){return arr[col][row];} //returns the (i,j)th entry
        bool operator<(const SO6 &) const;
        const Z2& operator()(int8_t i, int8_t j) const {return arr[i][j];} //returns the (i,j)th entry but for const
        bool operator==(SO6&); //checking equality up to signed permutation
        bool operator==(const SO6 &) const;
        Z2* operator[](const int8_t i) {return arr[i];}  // Return the array element needed. 
        const Z2* operator[](const int8_t i) const {return arr[i];}  // Return the array element needed. 
        SO6 tMultiply(const int &);
        SO6 tMultiplyTranspose(const int &);        
        SO6 tMultiply(const int &,const int &,const int &);
        SO6 pattern_mod();
        void genLDE(); //generates LDE, called after multiplication and constructor
        friend std::ostream& operator<<(std::ostream&,const SO6&); //display
        static const SO6 identity() {
            SO6 I;
            for(int8_t k =0; k<6; k++) {
                I(k,k) = 1;
            }
            return I;
        }
        int8_t getLDE();
        SO6 getPattern();
        SO6 transpose();
        std::string printName(); 
        void row_permute(int *);
        std::vector<int8_t> hist;
    private:
        Z2 arr[6][6];
};
