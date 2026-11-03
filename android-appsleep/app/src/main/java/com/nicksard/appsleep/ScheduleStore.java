package com.nicksard.appsleep;

import android.content.Context;
import android.content.SharedPreferences;

import org.json.JSONArray;
import org.json.JSONException;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Saves schedules as JSON in SharedPreferences and keeps a parsed copy in memory. */
public final class ScheduleStore {
    private static final String PREFS = "schedules";
    private static final String KEY = "list";

    private static List<Schedule> cache;

    private ScheduleStore() {}

    private static SharedPreferences prefs(Context c) {
        return c.getApplicationContext().getSharedPreferences(PREFS, Context.MODE_PRIVATE);
    }

    public static synchronized List<Schedule> load(Context c) {
        if (cache == null) {
            List<Schedule> list = new ArrayList<>();
            String raw = prefs(c).getString(KEY, "[]");
            try {
                JSONArray arr = new JSONArray(raw);
                for (int i = 0; i < arr.length(); i++) {
                    list.add(Schedule.fromJson(arr.getJSONObject(i)));
                }
            } catch (JSONException ignored) {
                // Corrupt data: start over rather than crash the blocker.
            }
            cache = Collections.unmodifiableList(list);
        }
        return cache;
    }

    public static synchronized Schedule find(Context c, String id) {
        if (id == null) return null;
        for (Schedule s : load(c)) if (s.id.equals(id)) return copy(s);
        return null;
    }

    public static synchronized void save(Context c, Schedule schedule) {
        List<Schedule> list = new ArrayList<>(load(c));
        boolean replaced = false;
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).id.equals(schedule.id)) {
                list.set(i, schedule);
                replaced = true;
            }
        }
        if (!replaced) list.add(schedule);
        write(c, list);
    }

    public static synchronized void delete(Context c, String id) {
        List<Schedule> list = new ArrayList<>(load(c));
        list.removeIf(s -> s.id.equals(id));
        write(c, list);
    }

    /** Returns the first schedule that has this app asleep right now, or null. */
    public static Schedule sleepingScheduleFor(Context c, String packageName, LocalDateTime now) {
        for (Schedule s : load(c)) {
            if (s.isAsleep(packageName, now)) return s;
        }
        return null;
    }

    private static Schedule copy(Schedule s) {
        try {
            return Schedule.fromJson(s.toJson());
        } catch (JSONException e) {
            throw new IllegalStateException(e);
        }
    }

    private static void write(Context c, List<Schedule> list) {
        JSONArray arr = new JSONArray();
        try {
            for (Schedule s : list) arr.put(s.toJson());
        } catch (JSONException e) {
            throw new IllegalStateException(e);
        }
        prefs(c).edit().putString(KEY, arr.toString()).apply();
        cache = Collections.unmodifiableList(list);
    }
}
