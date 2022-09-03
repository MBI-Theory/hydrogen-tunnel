.PHONY: check check-min check-dble check-quad clean

check:	check-dble

check-min:
	../hydrogen_tunnel_v2.x < h000-0_100-dble.inp > h000-0_100-dble.out
	./scripts/check.sh h000-0_100-dble.out
	../general_tunnel.x < gyu000-0_050-3c-out-dble.inp > gyu000-0_050-3c-out-dble.out
	./scripts/check.sh gyu000-0_050-3c-out-dble.out

check-dble:
	./scripts/run-all.sh dble

check-quad:
	./scripts/run-all.sh quad

clean:
	rm -f *.out *.table
