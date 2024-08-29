EFLAGS = -Wall -Wextra -Werror -std=c11
CFLAGS = $(EFLAGS) -c
OFLAGS = $(EFLAGS) -o
S21_MATRIX_C = s21_matrix.c
S21_MATRIX_O = s21_*.o

all: s21_matrix.a style test gcov_report

s21_matrix.a:
	gcc $(CFLAGS) $(S21_MATRIX_C)
	ar r s21_matrix.a $(S21_MATRIX_O)
	rm -rf *.o
	ranlib s21_matrix.a

style:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.c *.h tests/*.c
	rm .clang-format

test: s21_matrix.a
	gcc $(CFLAGS) tests/test_s21_matrix.c
	gcc $(OFLAGS) test test_s21_matrix.o s21_matrix.a -lcheck -lm
	./test

rebuild: clean test

gcov_report:
	gcc --coverage $(CFLAGS) $(S21_MATRIX_C)
	ar r s21_matrix.a $(S21_MATRIX_O)
	ranlib s21_matrix.a
	gcc $(CFLAGS) tests/test_s21_matrix.c
	gcc --coverage $(OFLAGS) test test_s21_matrix.o s21_matrix.a -lcheck -lm
	rm -rf *.o *.a
	./test
	gcovr -r . --html --html-details -o report.html
	open report.html

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt test.info test.dSYM proj