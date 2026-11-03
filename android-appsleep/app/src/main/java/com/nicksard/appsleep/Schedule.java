package com.nicksard.appsleep;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import java.time.DayOfWeek;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;
import java.time.format.FormatStyle;
import java.time.format.TextStyle;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.Set;
import java.util.UUID;

/** One sleep window: a set of apps that are asleep between two times on chosen days. */
public final class Schedule {
    public String id = UUID.randomUUID().toString();
    public String name = "Bedtime";
    /** Minutes after midnight. If start > end the window runs overnight. Equal means all day. */
    public int startMinute = 22 * 60;
    public int endMinute = 7 * 60;
    /** Bit (d - 1) is set when day d is on, where d is DayOfWeek.getValue() (Monday = 1). */
    public int daysMask = 0b1111111;
    public final Set<String> packages = new LinkedHashSet<>();
    public boolean enabled = true;
    public boolean silenceNotifications = true;

    public boolean isDayOn(DayOfWeek day) {
        return (daysMask & (1 << (day.getValue() - 1))) != 0;
    }

    public void setDayOn(DayOfWeek day, boolean on) {
        int bit = 1 << (day.getValue() - 1);
        daysMask = on ? (daysMask | bit) : (daysMask & ~bit);
    }

    /** True when this schedule is switched on and the given moment falls inside its window. */
    public boolean isActiveAt(LocalDateTime t) {
        if (!enabled || daysMask == 0) return false;
        int minute = t.getHour() * 60 + t.getMinute();
        DayOfWeek today = t.getDayOfWeek();
        if (startMinute == endMinute) {
            return isDayOn(today);
        }
        if (startMinute < endMinute) {
            return isDayOn(today) && minute >= startMinute && minute < endMinute;
        }
        // Overnight: the evening part belongs to today, the early-morning part to yesterday's window.
        return (isDayOn(today) && minute >= startMinute)
                || (isDayOn(today.minus(1)) && minute < endMinute);
    }

    public boolean isAsleep(String packageName, LocalDateTime t) {
        return packages.contains(packageName) && isActiveAt(t);
    }

    public static String formatMinute(int minuteOfDay) {
        LocalTime time = LocalTime.of(minuteOfDay / 60, minuteOfDay % 60);
        return time.format(DateTimeFormatter.ofLocalizedTime(FormatStyle.SHORT));
    }

    public String timeRangeLabel() {
        if (startMinute == endMinute) return "All day";
        return formatMinute(startMinute) + " – " + formatMinute(endMinute);
    }

    public String daysLabel() {
        if (daysMask == 0b1111111) return "Every day";
        if (daysMask == 0b0011111) return "Weekdays";
        if (daysMask == 0b1100000) return "Weekends";
        if (daysMask == 0) return "No days";
        StringBuilder sb = new StringBuilder();
        for (DayOfWeek d : DayOfWeek.values()) {
            if (isDayOn(d)) {
                if (sb.length() > 0) sb.append(", ");
                sb.append(d.getDisplayName(TextStyle.SHORT, Locale.getDefault()));
            }
        }
        return sb.toString();
    }

    public JSONObject toJson() throws JSONException {
        JSONObject o = new JSONObject();
        o.put("id", id);
        o.put("name", name);
        o.put("start", startMinute);
        o.put("end", endMinute);
        o.put("days", daysMask);
        o.put("enabled", enabled);
        o.put("silence", silenceNotifications);
        JSONArray apps = new JSONArray();
        for (String p : packages) apps.put(p);
        o.put("apps", apps);
        return o;
    }

    public static Schedule fromJson(JSONObject o) {
        Schedule s = new Schedule();
        s.id = o.optString("id", s.id);
        s.name = o.optString("name", s.name);
        s.startMinute = o.optInt("start", s.startMinute);
        s.endMinute = o.optInt("end", s.endMinute);
        s.daysMask = o.optInt("days", s.daysMask);
        s.enabled = o.optBoolean("enabled", true);
        s.silenceNotifications = o.optBoolean("silence", true);
        JSONArray apps = o.optJSONArray("apps");
        if (apps != null) {
            for (int i = 0; i < apps.length(); i++) s.packages.add(apps.optString(i));
        }
        return s;
    }
}
