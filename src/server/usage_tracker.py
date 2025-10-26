from __future__ import annotations
import os
import sys

IS_WEB_APP = os.environ.get("IS_WEB_APP")
USE_REDIS = bool(IS_WEB_APP)
REDIS_ERROR = "" if IS_WEB_APP else "Not running as WebApp"
if USE_REDIS:
    try:
        from redis import Redis
        from redis.exceptions import ConnectionError
        r = Redis.from_url(os.environ.get("REDIS_URL", ""))
        r.ping()
    except ImportError:
        USE_REDIS = False
        REDIS_ERROR = "Redis not installed"
        print(REDIS_ERROR, file=sys.stderr)
    except ValueError:
        USE_REDIS = False
        REDIS_ERROR = "No database set"
        print("No database set, set environment variable REDIS_URL to connect to a redis DB. Not needed for personal use.", file=sys.stderr)
    except ConnectionError:
        print("Can't connect to redis database, ping failed", file=sys.stderr) #fine to keep trying, may have been a one-off issue


def add_run_to_db(user_id: str, program_name: str):
    if USE_REDIS:
        try:
            #r.hincrby(f"user:{user_id}:programs", program_name, 1) #get runs per user
            r.hincrby(f"program:{program_name}:users", user_id, 1) #get users per program

            r.hincrby(f"program_runs", program_name, 1)  # get users per program

            r.sadd("all_users", user_id)  # unique user IDs
            r.sadd(f"program:{program_name}", user_id)  # users for each program
        except ConnectionError: #don't stop request if DB is down
            print("Can't connect to redis database for program result, connection failed.", file=sys.stderr)

def get_stats():
    if not USE_REDIS: return "Can't connect to database. " + REDIS_ERROR
    try:
        return f"""Total Users: {r.scard('all_users')}
TFOFinder: Users - {r.scard('program:TFOFinder')}, Runs - {int(r.hget('program_runs', 'TFOFinder') or 0)}
PinMol: Users - {r.scard('program:PinMol')}, Runs - {int(r.hget('program_runs', 'PinMol') or 0)}
smFISH: Users - {r.scard('program:smFISH')}, Runs - {int(r.hget('program_runs', 'smFISH') or 0)}"""
    except ConnectionError:
        print("Can't connect to redis database for stats, connection failed.", file=sys.stderr)
        return "Can't connect to database, connection failed"


